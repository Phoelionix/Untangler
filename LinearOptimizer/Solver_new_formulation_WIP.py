#===========================================================================
# Untangler: Free ensemble models from local minima with the wrong altlocs 
# Copyright (C)  2025 Spencer Passmore (spencerpassmore@swin.edu.au)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#===========================================================================

# Finds assignment of altlocs that minimizes the geometric "badness" defined in LinearOptimizer.Input. 
# It may be easier to think about it as connecting up the atoms in a way that minimizes the total energy.


# TODO 
# We will quickly reach an intractable number of variables as we consider more conformations.
# (for each angle, there will be n^3 connections.)
# So it will be important to break down the problem in some way. Possibly a stochastic way.
# Also need to look at exploring different branches of solutions in an intelligent way, and with parallel processing.  
# Meaning of constraint is different in Input.py, where it refers to constraints between (disordered) atoms, and Solver.py, where it refers to constraints between conformers (ordered atoms).

# CRUCIAL reduce weight of clashes and compensate with an elastic constraint for number of clashes



import pulp as pl
import numpy as np
from enum import Enum
from LinearOptimizer.Input import *
from LinearOptimizer.VariableID import *
import itertools
import UntangleFunctions
import json
from copy import deepcopy
import gc; 
import sys
import shutil
import matplotlib.pyplot as plt
import matplotlib
from typing import Union
from time import time,sleep
import random


ALTERNATIVE_CONSTRAINT_2_FORMULATION =True

#THREADS=None # Num cpu threads made available to ILP solver 
THREADS=UntangleFunctions.NUM_THREADS



#import pulp as pl
# just for residues for now
PLOTTING=True
MAX_BOND_CHANGES_SECOND_HALF_ONLY=False
CHANGES_MUST_INVOLVE=None#["A"] # In testing.
DEBUG_FIRST_100_SITES=False
FORBID_MULTIPLE_LOCAL_CHANGES_WHEN_MAIN_CHAIN_ONLY=True
ALLOW_ALL_POSITION_CHANGE_GEOMECTIONS= True # Only if modify_forbid_conditions = True
FORCE_ALT_COORDS=False

KEEP_PREVIOUS_FLEXI=False
NUM_RELEASE_ROUNDS=0

# Specify None to consider all altlocs.
#ALTLOC_RUN_SUBSET_SIZES=[3,3,3,3,4,4,4,5,5,None] # None 
#ALTLOC_RUN_SUBSET_SIZES=[2,2,2,3,3,3,3,4,4,4] # None 
#ALTLOC_RUN_SUBSET_SIZES=[2,2,2,2,3,3,4,4,4,4,4,5] # None 
#ALTLOC_RUN_SUBSET_SIZES=[4,4,5,6] # None 


if UntangleFunctions.NO_UNRESTRAINED:
    MIN_ALTLOCS_TO_GLUE_GOOD_GEOMETRY_GROUPS=6
    MIN_ALTLOCS_TO_FIX_RANDOM=6
else:
    MIN_ALTLOCS_TO_GLUE_GOOD_GEOMETRY_GROUPS=3
    MIN_ALTLOCS_TO_FIX_RANDOM=3

#ALTLOC_RUN_SUBSET_SIZES=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,None] # None 
#ALTLOC_RUN_SUBSET_SIZES=[None] # None 
ALTLOC_RUN_SUBSET_SIZES_AFTER_DIFFICULT=[3,3,3,4,5,None]

DIFFICULT_SOLVE_TIME_THRESHOLD_IN_MINS=99999




ALTLOC_RUN_SUBSET_SIZES=[3,3,3,4,4,4] # None 
MIN_ALTLOCS_TO_GLUE_GOOD_GEOMETRY_GROUPS=3
MIN_ALTLOCS_TO_FIX_RANDOM=3

# TODO cplex solution pool

###
BETTER_THAN_WORST=0
BETTER_THAN_AVERAGE=1
BETTER_THAN_AVERAGE_MULT=2
tolerate_score_mode=BETTER_THAN_WORST
#MIN_IMPROVEMENT_FACTOR_TO_TOLERATE=10  # Higher is stricter (fewer high sigma connections will be allowed)
###

MEMORYLIMITGB=35




def add_sos(lp_problem:LpProblem,sos_name,sos_rule):
    lp_problem.sos1[sos_name]=sos_rule
def add_sos2(lp_problem:LpProblem,sos2_name,sos2_rule):
    lp_problem.sos2[sos2_name]=sos2_rule


def solve(chunk_sites: list[AtomChunk],disordered_connections:dict[str,list[LP_Input.Geomection]],out_dir,out_handle:str,force_no_flips=False,num_solutions=20,force_sulfur_bridge_swap_solutions=False,
          inert_protein_sites=False,protein_sites:bool=True,water_sites:bool=True,max_mins_start=5,mins_extra_per_loop=0.1,#max_mins_start=100,mins_extra_per_loop=10,
          inert_water_sites=False,
          #gapRel=0.001,
          #gapRel=0,
          gapRel=0.005,
          #gapRel=0.03,
          forbid_altloc_changes={"name":[]}, forbidden_atom_bond_changes={"name":[]},forbidden_atom_any_connection_changes={"name":[]},
          MAIN_CHAIN_ONLY=False,SIDE_CHAIN_ONLY=False,NO_CB_CHANGES=False,NO_ISOLATED_O_BOND_CHANGES=False, # Forbids BOND changes that do not involve the specified atoms.
          #max_bond_changes=None):  
          #max_bond_changes=24):  
          #max_bond_changes=7):  
          max_bond_changes=None,
          modify_forbid_conditions=True,  # Whether the function is allowed to modify (generally, turn off) the "forbidden" flags of the connections
          change_punish_factor=0, # Adds cost of: change_punish_factor x num_conformer_labels_changed/num_conformers x original_cost. Num conformers excludes conformers that are being ignore (see `site_being_considered()`).
          forbid_ring_changes=False,
          forbid_solutions_composed_of_better_solutions=False,
          forbid_CECD12_changes=False, # Leaving true until fix bug from missing interchangability of C[E/D]1 and C[E/D]2 when reading geometry restraints 
          reference_pdb_file=None
          ):  
          #max_bond_changes=None):  
    print("****************\nConstructing ILP problem\n****************")

    # protein_sites, water_sites: Whether these can be swapped.
    # gaprel : relative gap tolerance for the solver to stop (fraction) # https://coin-or.github.io/pulp/technical/solvers.html
    assert MAIN_CHAIN_ONLY + SIDE_CHAIN_ONLY + NO_CB_CHANGES <=1 
    
    TOSYMBOL="~TO~" # >>
    CEMSYMBOL="~CEM~"    

    if gapRel == 0:
        gapRel = None

    connection_types= (
        ConstraintsHandler.BondConstraint,
        ConstraintsHandler.AngleConstraint,
        ConstraintsHandler.NonbondConstraint,
        ConstraintsHandler.ClashConstraint,#1e2,
        ConstraintsHandler.TwoAtomPenalty,
    )
    num_forbidden_connections={k:0 for k in connection_types}
    num_allowed_connections={k:0 for k in connection_types}
    num_allowed_alternative_connections={k:0 for k in connection_types}
    num_small_fry_geomections={k:0 for k in connection_types}

    # if inert_protein_sites:
    #     assert protein_sites
    # if inert_water_sites:
    #     assert water_sites
    #first_atom_id = "1.N"
    lowest_site_num = np.inf
    for chunk_site in chunk_sites:
        if chunk_site.get_site_num() < lowest_site_num:
            lowest_site_num = chunk_site.get_site_num()


    #nodes = [chunk.unique_id() for chunk in chunk_sites.values()]

    # def get_variables(unique_id)->dict:
    #     d={}
    #     d["depth"] = unique_id.split("&")[0]
    #     d["site_num"],d["altloc"] = unique_id.split("&")[1].split(".")
    #     d["name"] = unique_id.split("&")[2]
    #     return d
    
    # def get(unique_id,key):
    #     return get_variables(unique_id)[key]

    def site_being_considered(site:Union[VariableID,AtomChunk]):
        assert type(site) in [VariableID,AtomChunk]
        if DEBUG_FIRST_100_SITES and site.get_site_num()>100:
           return False
        return (((protein_sites and not site.is_water)
            or (water_sites and site.is_water ))
            )
            #and site.site_num < 10) #XXX DEBUG TEMPORARY   


    forced_swap_solutions=[]




    lp_problem = pl.LpProblem(f"Untangling_Problem-{out_handle}", LpMinimize)

    

   

    ##### Constraints ####



    #Disordered variable id
    # Messy...
    #bond_choices = {} # Each must sum to 1
    distance_vars = [] 
    # NOTE When assigning LP variable names, the "class" of variable should follow format of variableName_other_stuff   (class of variable is variable_type.split("_")[0]) 
    constraint_var_dict:dict[VariableID,tuple[LP_Input.Geomection,pl.LpVariable]] = {} 
    site_var_dict:dict[VariableID,dict[str,dict[str,pl.LpVariable]]] ={} # site label vars
    site_altposvar_dict:dict[VariableID,dict[str,dict[int,pl.LpVariable]]] ={} # site position change vars
    site_altpos_dict:dict[VariableID,dict[str,dict[int,Atom]]]={}
    var_dictionaries = dict(
        constraints = constraint_var_dict,
        sites = site_var_dict,
    )

    # Setup where atoms are assigned.
    # Track which are which based on original altloc.

    all_altlocs:list[str]=[]
    for ch in chunk_sites:
        if ch.altloc not in all_altlocs:
            all_altlocs.append(ch.altloc)
    # def site_variable_name(site: VariableID, from_altloc:str, to_altloc:str):
    #     return f"atomState_{site}_{from_altloc}.{to_altloc}"
    disordered_atom_sites:list[VariableID] = []

    # TODO replace this system with one that just adds the N best constraints in round 1, 2N best constraints in round 2, etc.
    # TODO longer term - specify priority order?
    #improvement_factors_to_tolerate=np.array([100,8,4,3,2,1.5,1.25,1,0.9,]) 
    # if tolerate_score_mode==BETTER_THAN_AVERAGE:
    #     if len(all_altlocs)==2:
    #         improvement_factors_to_tolerate=np.array([100,2,1]) 
    #     elif len(all_altlocs)<=6:
    #         improvement_factors_to_tolerate=np.array([100,3,2,1.5,1.25,1.1,1,0.95,0.9,0.85,0.8,0.75,0.7,0.65]) 
    #     else:
    #         improvement_factors_to_tolerate=np.array([100,6,3,2,1.5,1.25,1.175,1.1,1.05,1]) 
    #mainchain_improvement_factor_requirement_mult_bond=1#0.5
    #mainchain_improvement_factor_requirement_mult_angle=1#0.01 #0.01
    mainchain_improvement_factor_requirement_mult_bond=0.5
    mainchain_improvement_factor_requirement_mult_angle=0.01 #0.01
    sidechain_improvement_factor_requirement_mult_bond=1
    sidechain_improvement_factor_requirement_mult_angle=1
    #TODO try a round where everything costs zero except clashes?
    if True:
        if len(all_altlocs)==2:
            improvement_factors_to_tolerate=np.array([100,2,1]) 
        elif len(all_altlocs)<=4:
            #improvement_factors_to_tolerate=np.array([100,12,10,4,3,2,1.5,1.2,1.1,1,0.95]) 
            improvement_factors_to_tolerate=np.array([2,1]) 
        else:
            #improvement_factors_to_tolerate=np.array([100,12,10,8,6,5,4.5,4,3.5,3,2.5,2.25,2,1.75,1.5,1.35,1.2,1.1,1,0.95,0.90]) 
            #improvement_factors_to_tolerate=np.array([100,12,6,4,3,2,1.5,1,0.95,0.90,0.85]) 
            #improvement_factors_to_tolerate=np.array([5,1.75,1.2,1,0.90,0.85,0.80]) 
            #improvement_factors_to_tolerate=np.array([1.75,1.5,1.25,1,0.9,0.8]) 
            #improvement_factors_to_tolerate=np.array([1,0.8]) 
            #improvement_factors_to_tolerate=np.array([10,5,4,3,2,1.5,1]) 
            #improvement_factors_to_tolerate=np.array([1]) 
            #improvement_factors_to_tolerate=np.array([5,2,1,0.8,0.6,0.5]) 
            #improvement_factors_to_tolerate=np.array([5,3,2,1,0.8,0.6,0.5]) 
            if UntangleFunctions.NO_UNRESTRAINED:
                improvement_factors_to_tolerate=np.array([5,2,1,0.5,0.25]) 
            else:

                improvement_factors_to_tolerate=np.array([100,4,2,1.5,1,0.75]) 

    #TODO limit alternatives to consider to the top N alternatives. Otherwise when have really bad outliers, introduce a huge number of branches.
    # TODO dynamical solution space size. Stop solve if taking too long, and increase the required improvement_factor, then retry. 

    #improvement_factors_to_tolerate=np.array([0.7]) # TODO replace this system with one that just adds the N best constraints in round 1, 2N best constraints in round 2, etc.




    # Set up atom swap variables
    forbid_flip_site_altloc_dict:dict[VariableID,list[str]] = {}
    for i, chunk in enumerate(chunk_sites):
        site = VariableID.Atom(chunk)

        if not site_being_considered(site):
            continue

        if site not in disordered_atom_sites:
            disordered_atom_sites.append(site)

        #from_altloc = get(chunk.unique_id(),"altloc")
        from_altloc = chunk.altloc

        if chunk.echoed_altloc is not None:
            if site not in forbid_flip_site_altloc_dict:
                forbid_flip_site_altloc_dict[site]=[]
            forbid_flip_site_altloc_dict[site].append(from_altloc)

        if chunk.has_alternate_coords():
            if site not in site_altposvar_dict:
                site_altposvar_dict[site]={} # site position change vars
                site_altpos_dict[site]={} 
            site_altposvar_dict[site][from_altloc]={i:None for i in chunk.alt_pos_options}
            site_altpos_dict[site][from_altloc]=chunk.alt_pos_options
        
        #TODO!!!!!!! if forbidden due to being absent from any disordered connections that involve the site, exclude! 
        if site not in site_var_dict:
            site_var_dict[site] = {}
        site_var_dict[site][from_altloc]={}


    '''
    for site in disordered_atom_sites:
        if site.is_water:
            have_water=True
            break
    else:
        have_water=False

    def get_force_no_flips_name(site,altloc):
        return f"forceNoFlips_{site}_{altloc}"

    for site in disordered_atom_sites:
        if not site_being_considered(site) :
            continue
        site_altlocs = []
        for possible_altloc in site_var_dict[site]:
            site_altlocs.append(possible_altloc)


        allowed_to_altlocs=site_altlocs # NOTE May want to change this to be all altlocs if water is involved. 
        if len(all_altlocs)>2: #TODO optimize
            for from_altloc in site_altlocs:
                # Create variable for each possible swap to other altloc
                for to_altloc in allowed_to_altlocs:
                    var_atom_assignment =  pl.LpVariable(
                        f"@{site}_{from_altloc}.{to_altloc}",
                        lowBound=0,upBound=1,cat=pl.LpBinary #TODO pl.LpBinary
                    )
                    # For warm start.
                    var_atom_assignment.setInitialValue(0)
                    if to_altloc == from_altloc: 
                        var_atom_assignment.setInitialValue(1)
                    site_var_dict[site][from_altloc][to_altloc]=var_atom_assignment

                # Each ordered atom is assigned to one conformation from:to == n:1
                lp_problem += (  
                    #lpSum(site_var_dict[site][from_altloc])==1,
                    lpSum(site_var_dict[site][from_altloc].values())==1,
                    f"fromAltLoc_{site}.{from_altloc}"
                )

            # Each conformation is assigned no more than one ordered atom. 
            for to_altloc in allowed_to_altlocs: 
                to_altloc_vars:list[LpVariable] = []
                for from_alt_loc_dict in site_var_dict[site].values():
                    to_altloc_vars.append(from_alt_loc_dict[to_altloc])
                # from:to == 1:n
                if set(site_altlocs)==set(all_altlocs):
                    lp_problem += (  
                        lpSum(to_altloc_vars)==1,
                        f"toAltLoc_{site}{TOSYMBOL}{to_altloc}"
                    )
                # E.g. water sites which are not in all conformations.
                else:
                    lp_problem += (  
                        lpSum(to_altloc_vars)<=1,
                        f"toAltLoc_{site}{TOSYMBOL}{to_altloc}"
                    )

        else:
            var_flipped =  pl.LpVariable(
                f"@{site}_Flipped",
                lowBound=0,upBound=1,cat=pl.LpBinary #TODO pl.LpBinary
            )
            var_flipped.setInitialValue(0)
            var_not_flipped = 1-var_flipped 

            if len(site_altlocs)==1:  
                site_var_dict[site][site_altlocs[0]][site_altlocs[0]] = var_not_flipped
                other_altlocs = [a for a in all_altlocs if a != site_altlocs[0]]
                assert len(other_altlocs)==1, (site, "You may have passed only a single altloc/conformation to the optimizer")
                site_var_dict[site][site_altlocs[0]][other_altlocs[0]] = var_flipped
            else:
                site_var_dict[site][site_altlocs[0]][site_altlocs[0]] = var_not_flipped
                site_var_dict[site][site_altlocs[0]][site_altlocs[1]] = var_flipped
                site_var_dict[site][site_altlocs[1]][site_altlocs[0]] = var_flipped
                site_var_dict[site][site_altlocs[1]][site_altlocs[1]] = var_not_flipped
            # lp_problem += (  
            #     #lpSum(site_var_dict[site][from_altloc])==1,
            #     var_flipped+var_not_flipped==1,
            #     f"fromAltLoc_{site}_{from_altloc}"
            # )

        # Alternate position variables
        if site in site_altposvar_dict:
            alt_indices=[]
            for from_altloc in site_altposvar_dict[site]:
                for alt_idx in site_altposvar_dict[site][from_altloc]:
                    if alt_idx not in alt_indices:
                        alt_indices.append(alt_idx) # XXX stupid
            
            # Create variable for each set of position changes, and assign the variables to the from_altlocs that they involve
            for alt_idx in alt_indices:
                var_atom_pos =  pl.LpVariable(
                    f"altCoords_{site}.{alt_idx}", # When these variables are 1, they will correspond to the position changes that we need to make
                    lowBound=0,upBound=1,cat=pl.LpBinary 
                )
                var_atom_pos.setInitialValue(0)
                for from_altloc in site_altposvar_dict[site]:
                    if alt_idx in site_altposvar_dict[site][from_altloc]:
                        site_altposvar_dict[site][from_altloc][alt_idx]=var_atom_pos
            # 
            for from_altloc in site_altposvar_dict[site]:
                var_no_change =  pl.LpVariable(
                    f"sameCoords_{site}.{from_altloc}",
                    lowBound=0,upBound=1,cat=pl.LpBinary 
                )
                # Create "no change variable" that will be set to 1 when the alt position variables are 0 (due to following condition)
                var_no_change.setInitialValue(1) 
                site_altposvar_dict[site][from_altloc][0]=var_no_change
                # Only one position
                lp_problem += (  
                    lpSum(site_altposvar_dict[site][from_altloc].values())==1,
                    f"altCoords_{site}_{from_altloc}"
                )
                
                if FORCE_ALT_COORDS and modify_forbid_conditions:
                    lp_problem += (  
                        var_no_change==0,
                        f"forceAltCoords_{site}_{from_altloc}"
                    )

        need_anchor = ( 
            not (have_water and inert_water_sites)     # if water sites are inert they break the symmetry
            and change_punish_factor==0
        )

        
        no_flip_from_altlocs=[]
        if (force_no_flips
            or ((site.site_num == lowest_site_num) and need_anchor) # Anchor solution to one where first disordered atom is unchanged.  
            or (not site.is_water and inert_protein_sites)
            or (site.is_water and inert_water_sites)
            or (site.atom_name in forbid_altloc_changes["name"])
        ) :
            no_flip_from_altlocs=site_altlocs 
        elif site in forbid_flip_site_altloc_dict:
            no_flip_from_altlocs=forbid_flip_site_altloc_dict[site]
        for altloc in no_flip_from_altlocs:
            lp_problem += (
                site_var_dict[site][altloc][altloc]==1,
                get_force_no_flips_name(site,altloc)
            ) 
        
        if CHANGES_MUST_INVOLVE is not None:
            for a in site_altlocs:
                for b in site_altlocs:
                    if a == b: # always allow no change.
                        continue 
                    if (a not in CHANGES_MUST_INVOLVE) and (b not in CHANGES_MUST_INVOLVE):
                        lp_problem += (
                            site_var_dict[site][a][b]==0,
                            f"Forbid_{site}_{a}{TOSYMBOL}{b}"
                        )  
    '''

    # TODO improve terminology
    # A "connection" just refers to a group of atoms with a constraint assigned by LinearOptimizer.Input.
    # A disordered connection refers to all the *possible* groupings of these atoms.

    # The variables of one of flexible_allowed_constr_sets[i] and flexible_forbidden_constr_sets[i] will be active for each "tolerance round" 
    # Whether they are active in tolerance round `r` is given by flexible_allowed_constr_tranches[i][r]
    # Each index i corresponds to a particular geomection (one variable for each altloc the geomection could be assigned).
    flexible_allowed_constr_sets:list[list[tuple]]=[]
    flexible_forbidden_constr_sets:list[list[tuple]]=[]
    flexible_forbidden_constr_types:list[Type]=[]
    flexible_vars:list[LpVariable]=[]
    flexible_allowed_constr_tranches:list[list[bool]]=[] #  TODO better name
    flexible_bad_original_forbid_constrs:list[list]=[] # List of list of constraints forbidding an original geomection. Used if and only if original geomection is inactive and always_allow_original_tranches[i][r] = False
    flexible_bad_original_allowed_constrs:list[list]=[]
    flexible_bad_original_constr_types:list[Type]=[]
    always_allow_original_tranches:list[list[bool]]=[] #  
    flexible_bad_original_vars:list=[]

    initial_badness=0
    def add_constraints_from_disordered_connection(constraint_type:VariableKind,disordered_connection: list[LP_Input.Geomection],global_score_tolerate_threshold=0):
        # Rule: If all atom assignments corresponding to a connection are active,
        # then all those atoms must be swapped to the same assignment.
        nonlocal lp_problem  
        nonlocal initial_badness
        nonlocal flexible_allowed_constr_sets
        nonlocal flexible_forbidden_constr_sets
        nonlocal flexible_forbidden_constr_types
        nonlocal flexible_bad_original_forbid_constrs
        nonlocal flexible_bad_original_constr_types
        nonlocal always_allow_original_tranches
        nonlocal flexible_allowed_constr_tranches


        def forbid_change_conditions():
            if not modify_forbid_conditions:
                return False
            
            MCH=["N","CA","CB","C","O"]
            chunks = disordered_connection[0].atom_chunks
            #involves_main_chain=any(ch.name in MCH for ch in chunks)
            in_main_chain=all(ch.name in MCH for ch in chunks)

            
            if forbid_CECD12_changes:
                if constraint_type==VariableKind.Bond or constraint_type==VariableKind.Angle: 
                    for ch in chunks:
                        if ch.get_resname() in ["TYR","PHE"]:
                            if ch.name in ["CD1","CD2","CE1","CE2","CZ"]:
                                for conn in disordered_connection: # Fix the bug...
                                    conn.ts_distance=0
                                    conn.z_score=0
                                return True
            if constraint_type==VariableKind.Bond: 
                # TODO this allows water to change. But necessary.
                if MAIN_CHAIN_ONLY and not in_main_chain and (any(ch.get_resname() not in ["CYS","HOH"] for ch in chunks)):  # XXX tidy up and put in a separate python file for specifying what to optimize
                    return True
                if SIDE_CHAIN_ONLY and in_main_chain and (any(ch.get_resname() not in ["HOH",] for ch in chunks)):
                    return True
                if NO_CB_CHANGES and any(ch.name =="CB" for ch in chunks):
                    return True
                if forbid_ring_changes:
                    for ch in chunks:
                        if ch.get_resname() in ["TYR","PHE"]:
                            if ch.name in ["CD1","CD2","CE1","CE2","CZ"]:
                                return True
                        elif ch.get_resname()=="PRO":
                            if ch.name in ["CG"]: # Include CD?
                                return True
                        elif ch.get_resname()=="TRP":
                            if ch.name not in MCH\
                            and ch.name not in ["CG"]:
                                return True


            # Forbid changes that are costly to consider and don't seem to tangle
            #TESTING_DISABLE_CHANGES=["CD2","CE2","OH","NH2","NZ"]
            TESTING_DISABLE_CHANGES=[]
            for ch in disordered_connection[0].atom_chunks:
                if constraint_type==VariableKind.Bond:
                    #if ch.name in ["O","OH"]:
                    #if ch.name in ["O"]:

                    #if ch.name in ["O","CD2"]:
                    #if ch.name in ["O","CD2","NH2"] or ch.name[0]=="O":
                    #if ch.name in ["O","OH","OG","OG1","OD1","NZ"]:
                    #if ch.name[0]=="O":
                    # if NO_ISOLATED_O_BOND_CHANGES and ch.name in ["O"]:
                    #     return True
                    #TODO implement option turn off swaps of O or other "endpoint" atom to its ridden atom (C for O) when there are no other connection changes involving same C
                    
                    if ch.name in forbidden_atom_bond_changes["name"]:
                        return True
                        
                if (ch.name in forbidden_atom_any_connection_changes["name"]) or (ch.name in TESTING_DISABLE_CHANGES):
                    return True
            return False
        
        forbid_constraint_change=forbid_change_conditions()
        

        scores = [conn.ts_distance for conn in disordered_connection]
        # if max(scores)-min(scores) < required_cost_range_to_consider:
        #     num_small_fry_geomections+=1
        #     return
        min_score = min(scores)
        
        # If change is allowed and the difference is minor, don't bother optimizing for it
        small_fry = []
        
        absolute_small_fry_scale = True
        if absolute_small_fry_scale:
            # If difference in cost from lowest costing ordered connection of the disordered connection is tiny, don't bother optimizing for it. (small fry)
            # Small fry variables are not included in the cost function. 
            # TODO forbid next-best solutions from reusing same set of non-small fry connections. 
            #required_cost_range_to_consider=1.0e-2
            required_cost_range_to_consider=0
            assert required_cost_range_to_consider>=0
            small_fry_threshold = min_score+required_cost_range_to_consider
        else:
            small_fry_factor = 0.1
            assert small_fry_factor >=0
            small_fry_threshold=min_score*(1+small_fry_factor)
        small_fry = [conn for conn in disordered_connection if conn.ts_distance < small_fry_threshold]
        
        #small_fry = [conn for conn in disordered_connection if conn.ts_distance - min_score < required_cost_range_to_consider]
        nonlocal num_small_fry_geomections
        num_small_fry_geomections[disordered_connection[0].connection_type]+=len(small_fry)
        
            

        

        if len(small_fry)==len(disordered_connection) and not forbid_constraint_change:
            return
        




        site_altlocs_same=True
        sites = [VariableID.Atom(ch) for ch in disordered_connection[0].atom_chunks]
        for site in sites:
            if not site_being_considered(site):
                    return # Don't add this constraint

        # dicts indexed by code corresponding to from altlocs (e.g. "ACB" means connecting up site 1 altloc A, site 2 altloc C, site 3 altloc B)
        geomection_var_dict:dict[str,tuple[LP_Input.Geomection,LpVariable]]={} # indexed by from_altloc

        #assert site_altlocs_same, (disordered_connection[0].connection_type, len(connection_dict),n**m,n,m)

        PENALIZE_O_BOND_CHANGES=False
        if PENALIZE_O_BOND_CHANGES:
            # Increase weight of bonds with O by fraction of average of current altlocs
            if (constraint_type==VariableKind.Bond
                and "O" in [ch.name for ch in disordered_connection[0].atom_chunks]):
                    avg_current_distance = np.mean([conn.ts_distance for conn in disordered_connection if conn.original()])

                    frac=0.2
                    for conn in disordered_connection:
                        if not conn.original():
                            conn.ts_distance+=frac*avg_current_distance

        def get_tag(geomection:LP_Input.Geomection):
            from_ordered_atoms = ";".join([f"{ch.resnum}.{ch.name}.{ch.altloc}" for ch in geomection.atom_chunks])
            tag=from_ordered_atoms
            extra_tag=""
            if geomection.hydrogen_tag!="":
                extra_tag = "Htag_"+geomection.hydrogen_tag+"_"
            tag=extra_tag+tag+geomection.poschange_tag
            return tag
        


        for geomection in disordered_connection:
            if geomection.original():
                initial_badness+=geomection.ts_distance
                
        if modify_forbid_conditions:
            #TODO Moves this logic to Input.py
            original_scores=[conn.ts_distance for conn in disordered_connection if conn.original()]

            
            # Threshold below which alternatives will be considered if flagged as forbidden by Input.py.

            #local_score_tolerate_threshold=2*worst_no_change_score
            #local_score_tolerate_threshold=2*worst_no_change_score
            #local_score_tolerate_threshold=10*worst_no_change_score
            #local_score_tolerate_threshold=3*worst_no_change_score  # * len(altlocs)?
            #local_score_tolerate_threshold=2*worst_no_change_score  # * len(altlocs)?   # TODO This should be done in input when setting what is forbidfden.
            if len(original_scores)==0: # FIXME
                local_score_tolerate_threshold=0
            elif tolerate_score_mode==BETTER_THAN_WORST:
                local_score_tolerate_threshold=max(original_scores)  # * len(altlocs)?   # TODO This should be done in input when setting what is forbidfden.
            elif tolerate_score_mode==BETTER_THAN_AVERAGE:
                local_score_tolerate_threshold=np.mean(original_scores)  # * len(altlocs)?   # TODO This should be done in input when setting what is forbidfden.
            elif tolerate_score_mode==BETTER_THAN_AVERAGE_MULT:
                local_score_tolerate_threshold=min(np.mean(original_scores)*2,max(original_scores))
            always_tolerate_score_threshold = max(local_score_tolerate_threshold,global_score_tolerate_threshold) #worst_no_change_score*10+1e4

    
            always_tolerate_score_threshold_sequence=always_tolerate_score_threshold/improvement_factors_to_tolerate
            # If is sidechain
            #if any(name not in ["CB","C","CA","N","O"] for name in disordered_connection[0].atom_names):
            if any(name not in ["C","CA","N","O"] for name in disordered_connection[0].atom_names):
                if disordered_connection[0].connection_type is ConstraintsHandler.AngleConstraint:
                    always_tolerate_score_threshold_sequence/=sidechain_improvement_factor_requirement_mult_angle
                elif disordered_connection[0].connection_type is ConstraintsHandler.BondConstraint:
                    always_tolerate_score_threshold_sequence/=sidechain_improvement_factor_requirement_mult_bond
            else:
                if disordered_connection[0].connection_type is ConstraintsHandler.AngleConstraint:
                    always_tolerate_score_threshold_sequence/=mainchain_improvement_factor_requirement_mult_angle
                elif disordered_connection[0].connection_type is ConstraintsHandler.BondConstraint:
                    always_tolerate_score_threshold_sequence/=mainchain_improvement_factor_requirement_mult_bond
            ########
            tmp_scores= [conn.ts_distance for conn in disordered_connection]
            tmp_scores.sort()
            add_flexi_up_to=int(len(all_altlocs)*10)
            #Nth_best_threshold=tmp_scores[min(add_flexi_up_to-1,len(tmp_scores)-1)]
            Nth_best_threshold=np.inf
            del tmp_scores; del add_flexi_up_to
            ########
        
        # TODO if hydrogen tag, make copy parent.

        def get_altlocs_key(from_altlocs,position_option_indices=None):
            altlocs_key=''.join(from_altlocs)
            if position_option_indices is not None and not all(i==0 for i in position_option_indices):
                altlocs_key+=''.join([str(i) for i in position_option_indices])
            return altlocs_key
        for geomection in disordered_connection:
            tag = get_tag(geomection)
            ########

            # if not allowed: 
            #     for to_altloc in all_altlocs:
            #         assignment_vars = [site_var_dict[VariableID.Atom(ch)][ch.get_altloc()][to_altloc] for ch in geomection.atom_chunks]
            #         lp_problem += (
            #             lpSum(assignment_vars) <=  len(assignment_vars)-1,   
            #             f"FORBID{constraint_type.value}_{tag}>>{to_altloc}"
            #             )
            #else:

            var_active = pl.LpVariable(f"{constraint_type.value}_{tag}",  #TODO cat=pl.LpBinary
                                lowBound=0,upBound=1,cat=pl.LpBinary)
            
            
            var_active.setInitialValue(0) 
            if geomection.original():
                var_active.setInitialValue(1)
            constraint_var_dict[VariableID(tag,constraint_type.value)]=(geomection,var_active)
            #group_vars.append(var_active)
            altlocs_key=get_altlocs_key(geomection.from_altlocs,
                                        geomection.position_option_indices if geomection.involves_position_changes() else None)

            assert altlocs_key not in geomection_var_dict, (geomection.get_disordered_connection_id(),altlocs_key,list(geomection_var_dict.keys()))
            geomection_var_dict[altlocs_key]=(geomection,var_active)
        del tag



        nonlocal num_allowed_connections
        nonlocal num_forbidden_connections

        num_original_connections = len([conn for conn,_ in geomection_var_dict.values() if conn.original()])

        # TODO refactor, make geomection_var_dict same as allowed_geomection_var_dict. Deal with allowed and not allowed in separate loops.

        allowed_geomection_var_dict:dict[str,tuple[LP_Input.Geomection,LpVariable]]={} # indexed by from_altloc
        for altlocs_key, (geomection, var_active) in geomection_var_dict.items():

            is_flexible=False
            is_bad_original_constraint=False
            if modify_forbid_conditions:
                allowed_sequence=[geomection.ts_distance<=threshold for threshold in always_tolerate_score_threshold_sequence]
                allowed_sequence=[(geomection.ts_distance<=Nth_best_threshold and allowed) for allowed in allowed_sequence]
                if geomection.original():
                    #  So possible solution is guaranteed, always allow connections of original structure.  
                    if geomection.forbidden and not all(allowed_sequence):
                        is_bad_original_constraint=True
                        always_allow_original_tranches.append(allowed_sequence)
                    else:
                        allowed=True 
                elif forbid_constraint_change:
                    allowed=False
                elif not geomection.forbidden:
                    allowed=True
                else:
                    if all(allowed_sequence):
                        allowed=True
                    elif not any(allowed_sequence):
                        allowed=False
                    else:
                        is_flexible=True
                        flexible_allowed_constr_tranches.append(allowed_sequence)
                    
                # chance_allow_anyway=0
                # if not allowed and chance_allow_anyway>0:
                #     if np.random.rand()<chance_allow_anyway:
                #         allowed=True

                
                if geomection.involves_position_changes() and ALLOW_ALL_POSITION_CHANGE_GEOMECTIONS:
                    allowed=True
            else:
                allowed = not geomection.forbidden
            
            assert is_flexible + is_bad_original_constraint <=1
            
                
            tag = get_tag(geomection)

                
            # Geomections contain ordered atoms from any (and likely multiple) from_altloc labels that 
            # *are to be assigned to the same conformation (to_altloc label)*
            assignment_vars:list[LpVariable]=[]
            if len(geomection.atom_chunks)==2:                   
                    pass
            else:
                # Angle geomection is active if 2 specific bond geomections are active.
                assert geomection.connection_type in [ConstraintsHandler.AngleConstraint,] # Can add other 3+ site constraints later
                
                def get_altlocs_key_or_echoed_altlocs_key(from_altlocs,position_option_indices,echoed_altlocs,d_id):
                    bond_altlocs_key=get_altlocs_key(
                        from_altlocs,
                        position_option_indices,
                    )
                    if bond_altlocs_key not in ALL_mega_geomection_var_dict[d_id]:
                        assert False
                        # Bond is between echoed sites, so get the bond constraint between the parent altlocs that the sites are echoes of.
                        assert not all(ea is None for ea in echoed_altlocs)
                        parent_bond_from_altlocs=[]
                        for echoed_altloc, from_altloc in zip(echoed_altlocs,from_altlocs):
                            parent_bond_from_altlocs.append(echoed_altloc if echoed_altloc is not None else from_altloc)
                        bond_altlocs_key=get_altlocs_key(
                            parent_bond_from_altlocs,
                            position_option_indices,
                        )
                        assert bond_altlocs_key in ALL_mega_geomection_var_dict[d_id], (d_id,bond_altlocs_key,from_altlocs,parent_bond_from_altlocs,echoed_altlocs,position_option_indices)

                    return bond_altlocs_key
                
                bond_vars=[]
                for i,(ch_A,ch_B) in enumerate(zip(geomection.atom_chunks[:-1],geomection.atom_chunks[1:])):
                    site_tags=(ch_A.get_disordered_tag(),ch_B.get_disordered_tag())
                    d_id = LP_Input.Geomection.construct_disordered_connection_id(
                        ConstraintsHandler.BondConstraint,site_tags
                    )                    
                    ordered_conn_from_altlocs = geomection.from_altlocs[i:i+2]
                    position_option_indices=geomection.position_option_indices[i:i+2]
                    echoed_altlocs = [ch.echoed_altloc for ch in geomection.atom_chunks[i:i+2]] 
                    if d_id not in ALL_mega_geomection_var_dict:
                        d_id = LP_Input.Geomection.construct_disordered_connection_id(
                        ConstraintsHandler.BondConstraint,list(reversed(site_tags))
                        )
                        ordered_conn_from_altlocs=list(reversed(ordered_conn_from_altlocs))
                        position_option_indices = list(reversed(position_option_indices))
                        echoed_altlocs=list(reversed(echoed_altlocs))
                    
                    try:
                        bond_altlocs_key=get_altlocs_key_or_echoed_altlocs_key(
                            ordered_conn_from_altlocs,
                            position_option_indices,
                            echoed_altlocs,
                            d_id
                        )
                        bond_vars.append(ALL_mega_geomection_var_dict[d_id][bond_altlocs_key][1])
                    except Exception as e:
                        print(e)
                        print(ordered_conn_from_altlocs)
                        print(bond_altlocs_key, list(ALL_mega_geomection_var_dict[d_id].keys()))
                        assert False
                assignment_vars=bond_vars

                    
            # if variable is inactive, cannot have all atoms assigned to the same altloc.
            # Note that since every connection option is looped through, this also means 
            # that if variable is active, all atoms will be assigned to the same altloc.
            
            
            if is_flexible or is_bad_original_constraint or allowed:
                allowed_geomection_var_dict[altlocs_key]=(geomection, var_active)
                if geomection not in small_fry:
                    assert geomection.ts_distance>=0, (geomection,geomection.ts_distance)
                    #if geomection.connection_type == ConstraintsHandler.ClashConstraint:
                    distance_vars.append(geomection.ts_distance*var_active)
            # TODO change "not (is_flexible or is_bad_original_constraint or allowed)" to "always_forbidden" variable.
            if not is_flexible and not is_bad_original_constraint:
                num_allowed_connections[geomection.connection_type]+=allowed
                if not geomection.original(): 
                    num_allowed_alternative_connections[geomection.connection_type]+=allowed 
                num_forbidden_connections[geomection.connection_type]+=not allowed 





            #### CONSTRAINT 1 "If all atoms are active in a connection, that connection is active" #####
            flexible_allowed_constrs=[]
            flexible_forbidden_constrs=[]
            flexible_original_allowed=[]
            flexible_original_forbidden=[]
           
            #assert len(geomection.atom_chunks)==len(assignment_vars) or geomection.involves_position_changes()   
            allowed_constraint=None
            if len(assignment_vars)>0:
                allowed_constraint = (
                    lpSum(assignment_vars) <=  len(assignment_vars)-1+var_active,   # Active if all assignment vars active.
                    f"ALLOW{constraint_type.value}_{tag}"
                )                      
            forbidden_constraint = (
                lpSum(assignment_vars) <=  len(assignment_vars)-1,   
                f"FORBID{constraint_type.value}_{tag}"
            )  
            if is_flexible:
                flexible_allowed_constr_sets.append(allowed_constraint)
                flexible_forbidden_constr_sets.append(forbidden_constraint)
                flexible_forbidden_constr_types.append(geomection.connection_type)
                flexible_vars.append(var_active)
            elif is_bad_original_constraint:
                flexible_bad_original_allowed_constrs.append(allowed_constraint)
                flexible_bad_original_forbid_constrs.append(forbidden_constraint)
                flexible_bad_original_constr_types.append(geomection.connection_type)
                flexible_bad_original_vars.append(var_active)
            else: # Whether allowed or not will depend on the round within a solution loop, as governed by `flexible_allowed_constr_tranches`  
                if allowed:
                    if allowed_constraint is not None:
                        lp_problem += allowed_constraint
                else: 
                    lp_problem += forbidden_constraint
                del allowed

        
        #geomection_var_dict = allowed_geomection_var_dict
        ## CONSTRAINT 2 "Num ordered geometries (i.e. 'connections') per disordered geometry must be unchanged"
        use_constraint_2=(geomection.connection_type is ConstraintsHandler.BondConstraint  # This might be dodgy
            and not ALTERNATIVE_CONSTRAINT_2_FORMULATION )
        ##### TODO  assess if need this
        enforce_angle_count_conservation=True
        if ((geomection.connection_type is ConstraintsHandler.AngleConstraint) 
            and enforce_angle_count_conservation):
            use_constraint_2=True
        ####################
        if use_constraint_2:
            dID = disordered_connection[0].get_disordered_connection_id()
            lp_problem += (
                lpSum([var_active for (_,var_active) in allowed_geomection_var_dict.values()])==num_original_connections,  # less than number of FROM altlocs. i.e. number of conformations it's currently involved in. 
                f"{dID}_{num_original_connections}_connections"
            )

        # # Debugging
        # var_track_num_active = pl.LpVariable(f"Debug_{constraint_type.value}_{tag}_connections",  #TODO cat=pl.LpBinary
        #             lowBound=0,upBound=99,cat=pl.LpInteger)
        # var_track_num_active.setInitialValue(num_conformations_involved_in_preswap)
        # lp_problem += (
        #     var_track_num_active ==  lpSum([var_active for (_,var_active) in geomection_var_dict.values()]),   
        #     f"Tracking_{constraint_type.value}_{tag}_connections"
        # )
        # distance_vars.append(1e3*var_track_num_active)
            

        return allowed_geomection_var_dict,geomection_var_dict
        
    worst_connection_before_swap=None
    worst_global_no_change_score=0
    for connection_id, ordered_connection_choices in disordered_connections.items():
        for c in ordered_connection_choices:
            if c.original() and c.ts_distance> worst_global_no_change_score:
                worst_global_no_change_score = c.ts_distance
                worst_connection_before_swap = c

        
    # TODO make it the 25th percentile or something.
    # TODO make it adapt based on solver time.
    global_score_tolerate_threshold=worst_global_no_change_score/100
    global_score_tolerate_threshold=0
    print("global score tolerate threshold:",global_score_tolerate_threshold)
    #global_score_tolerate_threshold=0
    
    #TODO this should replace 'constraint_var_dict'
    mega_geomection_var_dict:dict[str,dict[str,tuple[LP_Input.Geomection,LpVariable]]]={}
    ALL_mega_geomection_var_dict:dict[str,dict[str,tuple[LP_Input.Geomection,LpVariable]]]={}
    for i, (connection_id, ordered_connection_choices) in enumerate(disordered_connections.items()):
        if i % 250 == 0:
            print(f"Adding constraints {i}/{len(disordered_connections)} ({connection_id})")
        constraint_type = VariableKind[connection_id.split('_')[0]]  #XXX ?????
        result = add_constraints_from_disordered_connection(constraint_type,ordered_connection_choices,global_score_tolerate_threshold=global_score_tolerate_threshold)
        if result is not None:
            #XXX
            disordered_geomection_var_dict,disordered_geomection_var_dict_including_forbidden=result
            mega_geomection_var_dict[connection_id]=disordered_geomection_var_dict
            ALL_mega_geomection_var_dict[connection_id]=disordered_geomection_var_dict_including_forbidden
    print(f"Num allowed geomections: {num_allowed_connections}")
    print(f"Num allowed alternative geomections: {num_allowed_alternative_connections}")
    print(f"Num forbidden geomections: {num_forbidden_connections}")
    print(f"Num small fry: {num_small_fry_geomections}")

    if ALTERNATIVE_CONSTRAINT_2_FORMULATION:
        def constr2_name(var,other_var):
            return f"mutualExclu_{var.name}@{other_var.name}"
        for itr, disordered_geomection_var_dict in enumerate(mega_geomection_var_dict.values()):
            if itr%250==0:
                print(f"{itr}/{len(mega_geomection_var_dict)}")
            for constraint,var in disordered_geomection_var_dict.values():
                #if constraint.connection_type!=ConstraintsHandler.BondConstraint: # Might be better.
                if constraint.connection_type==ConstraintsHandler.AngleConstraint:
                    break
                for other_constraint, other_var in disordered_geomection_var_dict.values():
                    if other_constraint==constraint:
                        continue
                    if constr2_name(other_var,var) in lp_problem.constraints:
                        continue
                    not_shareable=False
                    for i in range(len(constraint.from_altlocs)):
                        if (constraint.from_altlocs[i]==other_constraint.from_altlocs[i]):
                            #constraint.position_option_indices[i]==other_constraint.position_option_indices[i]):
                            not_shareable=True
                            break
                    if not_shareable:
                        lp_problem+=(
                            var+other_var<=1,
                            constr2_name(var,other_var)
                        )

    max_bond_changes_tuple=None
    exclude_alt_positions_from_max_changes=True
    if max_bond_changes is not None or MAIN_CHAIN_ONLY:
        
        sys.setrecursionlimit(int(1e4)) 
        no_change_vars=[] # List of variables for whether no changes in a particular disordered connection
        no_change_vars_dict:dict[str,LpVariable]={}
        for connection_id, geomection_var_dict in mega_geomection_var_dict.items():
            constraint_type = VariableKind[connection_id.split('_')[0]]
            if constraint_type!=VariableKind.Bond:
                continue
            current_vars=[]
            excluded_vars=[]
            for conn,var in geomection_var_dict.values():
                if conn.original():
                    current_vars.append(var)
                if (exclude_alt_positions_from_max_changes and conn.involves_position_changes()):
                    excluded_vars.append(var)

                
            no_change_var = pl.LpVariable(f"NoChanges_{connection_id}",
                lowBound=0,upBound=1,cat=pl.LpBinary)
            no_change_var.setInitialValue(1)

            # lpSum(current_vars)+lpSum(excluded_vars) is the number that count as active (think of excluded_vars as "wildcards")
            # len(current_vars) is the number we need active.
            lp_problem+=(no_change_var*len(current_vars)<= lpSum(current_vars)+lpSum(excluded_vars),  
                        f"NoChangesRestraintInactive_{connection_id}")
            
            lp_problem+=(no_change_var>= lpSum(current_vars)+lpSum(excluded_vars)-len(current_vars)+1,  
                        f"NoChangesRestraintActive_{connection_id}")
            

                #f"NoChanges_{connection_id}"
        
            #lp_problem += no_change_var
            #no_change_vars.append(no_change_var)
            no_change_vars.append(no_change_var)
            no_change_vars_dict[connection_id] = no_change_var
    if max_bond_changes is not None:    
        max_bond_changes_tuple = (
            lpSum(no_change_vars)>=(len(no_change_vars)-max_bond_changes),
            f"Max{max_bond_changes}BondChanges"
        )
        def limit_bond_changes():
            nonlocal lp_problem
            nonlocal max_bond_changes_tuple
            print(f"limiting num bond changes to {max_bond_changes}")
            lp_problem += max_bond_changes_tuple
            max_bond_changes_tuple=None
        if not MAX_BOND_CHANGES_SECOND_HALF_ONLY:
            limit_bond_changes()

    if NO_ISOLATED_O_BOND_CHANGES:
        # Forbid changes to angles where the only change is through an oxygen bond change

 
        for connection_id, geomection_var_dict in mega_geomection_var_dict.items():
            constraint_type = VariableKind[connection_id.split('_')[0]]
            if constraint_type!=VariableKind.Bond:
                continue
            reference_conn = list(geomection_var_dict.values())[0][0]
            resnum = reference_conn.res_nums[0] # Anchor to first atom
            # if not "C" in [ach.name for ach in reference_conn.atom_chunks]:
            
            c_beta=False
            if not all([ach.name in ["O","C"] for ach in reference_conn.atom_chunks]):
                if all([ach.name in ["OG","CB"] for ach in reference_conn.atom_chunks]):
                    c_beta=True
                else:
                    continue    
            altpos_vars=[]
            
            isolated_altpos_changes_okay=False
            if isolated_altpos_changes_okay:
                for ch in reference_conn.atom_chunks:
                    site = VariableID.Atom(ch)
                    if site not in site_altposvar_dict:
                        continue
                    altpos_vars.extend([v for k, v in site_altposvar_dict[site][ch.get_altloc()].items() if k != 0])
            for altlocs_key, (CO_bond, CO_bond_var) in geomection_var_dict.items():
                if not CO_bond.original():
                    continue
                from_altloc = CO_bond.from_altlocs[0]
                
                # Get all other bonds involving the C that are same altloc
                other_bond_nochange_vars=[]
                other_bond_list = [(("N","C"),(resnum+1,resnum)), (("CA","C"),(resnum,resnum)),(("C","CA"),(resnum,resnum)), (("C","N"),(resnum,resnum+1))]
                if c_beta:
                    other_bond_list = [(("CA","CB"),(resnum,resnum))]
                for atom_names,resnums in other_bond_list: # TODO shouldn't need to do this. Should be sorted in Geomection class
                    disordered_tags= [DisorderedTag(num,name) for num,name in zip(resnums,atom_names)]
                    other_d_id= LP_Input.Geomection.construct_disordered_connection_id(ConstraintsHandler.BondConstraint,disordered_tags)
                    if other_d_id in mega_geomection_var_dict and altlocs_key in mega_geomection_var_dict[other_d_id]: # XXX
                        _other_bond, other_var = mega_geomection_var_dict[other_d_id][altlocs_key]
                        other_bond_nochange_vars.append(other_var)
                assert len(other_bond_nochange_vars)<=2, (len(other_bond_nochange_vars),other_bond_nochange_vars,other_d_id)
                bond_name = "CO"
                if c_beta:
                    bond_name = "CBOG"
                lp_problem +=  (
                    CO_bond_var>=1+lpSum(other_bond_nochange_vars)-len(other_bond_nochange_vars) - lpSum(altpos_vars),
                    f"forbidIsolated{bond_name}BondChangeRes{resnum}.{from_altloc}"
                )

        #print(resnum_no_change_vars_dict[4])

  


    if MAIN_CHAIN_ONLY and FORBID_MULTIPLE_LOCAL_CHANGES_WHEN_MAIN_CHAIN_ONLY:
        # Forbid multiple changes in bonds that connect separate residues (N, C, CA, and CB/SG in CYS) within residue and its neighbours (Linked CYS does not count as neighbour as of writing)

        # Point of focusing main chain is to get long range traps out.
        # Thus we want to forbid very short range differences. 
        max_resnum=-1
        resnum_no_change_vars_dict:dict[int,list[LpVariable]]={}
        for connection_id, geomection_var_dict in mega_geomection_var_dict.items():
            constraint_type = VariableKind[connection_id.split('_')[0]]
            if constraint_type!=VariableKind.Bond:
                continue
            reference_conn = list(geomection_var_dict.values())[0][0]
            resnum = reference_conn.res_nums[0] # Anchor to first atom
            #if not all(ach.name in ["N","CA","C","CB","SG"] for ach in reference_conn.atom_chunks):
            if not any(ach.name in ["N","CA","C","CB"] for ach in reference_conn.atom_chunks):
                continue
            if resnum not in resnum_no_change_vars_dict:
                resnum_no_change_vars_dict[resnum]=[]
                max_resnum = max(resnum,max_resnum)
            resnum_no_change_vars_dict[resnum].append(no_change_vars_dict[connection_id])
        #print(resnum_no_change_vars_dict[4])

        for resnum in range(1,max_resnum+1):
            # consider vars for bonds between N, C, CA in residue and neighbouring residues
            no_change_vars=[v for v in resnum_no_change_vars_dict[resnum]]
            if resnum > 1:
                no_change_vars.extend(resnum_no_change_vars_dict[resnum-1])
            if resnum < max_resnum:
                no_change_vars.extend(resnum_no_change_vars_dict[resnum+1])
            max_changes=1
            lp_problem +=  (
                lpSum(no_change_vars)>=(len(no_change_vars)-max_changes),
                f"Max{max_changes}BondChangesNearRes{resnum}"
            )

    # for disordered_geomection_var_dict in mega_geomection_var_dict.values():
    #     active_constraints=[(constraint,var) for constraint,var in disordered_geomection_var_dict.values() if var.value()>0]
    #     if len(active_constraints)==0:
    #         print("warning: , initial constraint inactive")
    #         continue
    #     if active_constraints[0][0].connection_type!=ConstraintsHandler.BondConstraint:
    #         continue
    #     for site in active_constraints[0][0].atom_chunks:
    #         if site.name not in ["N","CA","CB","C","O"]:
    #             main_chain_atoms=False
    #             break
    #     else:
    #         main_chain_atoms = True

    #     if not main_chain_atoms:
    #         continue
        
    #     no_change_var = (
    #         lpSum(constr[1] for constr in active_constraints)==len(active_constraints)
    #     )
    #     #lp_problem += no_change_var
    #     #no_change_vars.append(no_change_var)
    #     no_change_vars.append(no_change_var)
    
    # print(len(no_change_vars))
    # print(no_change_vars[0])
    # print(no_change_vars[-1])
    # lp_problem += (
    #     lpSum(no_change_vars)>=(len(no_change_vars)-max_main_chain_bond_changes),
    #     f"Max{max_main_chain_bond_changes}BondChanges"
    # )


    # if  force_sulfur_bridge_swap_solutions \
    #     and [ch.element for ch in connection.atom_chunks]==["S","S"]:
    #     forced_swap_solutions.append(
    #         var_dict[atom_a]["flipped"]==var_dict[atom_b]["flipped"]==1
    #     )

    if  force_sulfur_bridge_swap_solutions:
        raise Exception("force swap solutions unimplemented")



    if change_punish_factor!=0:
        raise NotImplementedError()
        num_conformers=len([site_being_considered(chunk) for chunk in chunk_sites])
        original_cost = lpSum([dist for dist in distance_vars]).value()

        
        change_punish_factor/num_conformers * original_cost
        for var in conformer_change_vars: # TODO
            distance_vars.append(change_punish_factor/num_conformers * original_cost*var)



    badness = lpSum([dist for dist in distance_vars])
    lp_problem += (
    badness,
    "badness",)
    #initial_badness = badness.value() 
    print(f"Initial badness: {initial_badness}")
        
    
    # TODO: Suspect we want the two badness to be close to equal. Because if they aren't similar badness, we wouldn't expect to see both? Unless some sort of structural/conformational "switch" in tension?

    
    ######################
    log_out_dir = os.path.join(out_dir,"linear_optimizer_logs",out_handle,"")
    if os.path.isdir(log_out_dir):
        shutil.rmtree(log_out_dir,ignore_errors=True)
    os.makedirs(log_out_dir)

    log_file = f"{log_out_dir}/HighLevelLog.Log"
    if not os.path.exists(log_file):
        with open(log_file,'w') as f:
            f.write(f"{out_handle} altloc optimizer log\n")
    def log(line:str):
        with open(log_file, 'a') as f:
            f.write(str(line)+"\n")
            

    swaps_file =  swaps_file_path(out_dir,out_handle,all_altlocs)
    if os.path.isfile(swaps_file):
        shutil.move(swaps_file,swaps_file+"#")
    def update_swaps_file(distances, site_assignment_arrays,record_notable_improvements_threshold=None): # record_notable_improvements_threshold: fractional improvement required to record separately
        assert len(distances)==len(site_assignment_arrays)
    # Create json file that lists all the site *changes* that are required to meet the solution. 
        out_dict = {"target": out_handle,"initial badness":initial_badness,"solutions":{}}
        if record_notable_improvements_threshold is not None:
            assert 0 <= record_notable_improvements_threshold < 1
            sep_dict= deepcopy(out_dict)
            best_improvement = 0
        else: 
            sep_dict=None
        for i, (distance, atom_assignments) in enumerate(zip(distances,site_assignment_arrays)):
            solution_dict = {"badness": distance}
            out_dict["solutions"][f"solution {i+1}"] = solution_dict
            if record_notable_improvements_threshold is not None and 1-distance/initial_badness >= record_notable_improvements_threshold: 
                best_improvement=max(best_improvement,1-distance/initial_badness)
                sep_dict["solutions"][f"solution {i+1}"]=solution_dict
            moves:dict[str,dict[str]] = {}
            solution_dict["moves"]=moves
            for site in atom_assignments:
                site_key = f"site {site}"
                for from_altloc,to_altloc in atom_assignments[site].items():
                    if from_altloc == to_altloc:
                        continue # No change needed
                    if site_key not in moves:
                        moves[site_key] = {}
                    moves[site_key][from_altloc] = to_altloc
        with open(swaps_file,'w') as f: 
            json.dump(out_dict,f,indent=4)


        if sep_dict is not None and len(sep_dict["solutions"])>0:
            separate_record_file = f"{log_out_dir}/Diff_{f'{best_improvement*100:.2f}'}.json"
            with open(separate_record_file,'w') as f2: 
                json.dump(sep_dict,f2,indent=4)


    assert num_solutions >= len(forced_swap_solutions) 

    site_assignment_arrays:list[dict[VariableID,dict[str,str]]]=[]
    distances=[]

    def write_current_connections(out_file):
        connections_str = ""
        connections_str+="name, ideal,model, z-score, cost\n"
        vals = [(constraint,var) for constraint,var in constraint_var_dict.values() if var.value()>0.5]
        #vals.sort(key=lambda x: x[0].z_score,reverse=True)
        vals.sort(key=lambda x: x[0].ts_distance,reverse=True)
        vals_selection = []
        ignore_zero_constraint_types = [ConstraintsHandler.ClashConstraint,ConstraintsHandler.NonbondConstraint,ConstraintsHandler.TwoAtomPenalty] 
        for constraint,var in vals:
            if constraint.ts_distance != 0 or (constraint.connection_type not in ignore_zero_constraint_types):
                vals_selection.append((constraint,var))
        vals = vals_selection

        for constraint, var in vals:
            connections_str+=f"{var.name} {constraint.ideal:.2e} {constraint.actual:.2e} {constraint.z_score:.2e} {constraint.ts_distance:.2e}\n"
        with open(out_file,'w') as f:
            f.write(connections_str)
    write_current_connections(f"{log_out_dir}/OriginalConnections.txt")

    create_initial_variable_files=True
    bonds_replaced_each_loop:list[list[LP_Input.Geomection]]=[]

    if create_initial_variable_files:
        #TODO make function with file name as arg.
        with open(f"{log_out_dir}/ProblemStatusBeforeFirstLoop.txt",'w') as f:
            f.write(f"Status: {LpStatus[lp_problem.status]}\n")

            for v in lp_problem.variables():
                try:
                    if v.value() > 0.5:
                        f.write(f"{v.name} = {v.value()}\n")
                except:
                    raise Exception(v,v.value())
            f.write(f"Total distance = {value(lp_problem.objective)}")


    inverted_active_variables = []
    inverted_active_variable_names = []
    non_original_variable_history = []
    previously_used_variable_names=[]

    flexifix_constr_names=[]
    def get_num_tranches():
        return len(improvement_factors_to_tolerate)
    def get_num_rounds():
        if get_num_tranches()==1:
            return 1
        return get_num_tranches() + KEEP_PREVIOUS_FLEXI*NUM_RELEASE_ROUNDS
    
    def flexifix_variables(r:int):
        # TODO could consider making this an elastic constraint.
        nonlocal flexifix_constr_names
        nonlocal lp_problem
        if r < get_num_tranches():
            for var in flexible_vars:
                flexifix_name=f"flexiFixed_{var.name}"
                if flexifix_name in flexifix_constr_names:
                    continue
                if value(var) > 0.5:
                    flexifix_constr_names.append(flexifix_name)
                    lp_problem+=(
                        var==1,
                        flexifix_name
                    )
        else: # unfix everything on final rounds
            print("Releasing flexi variables")
            if NUM_RELEASE_ROUNDS==0:
                warnings.warn("This shouldn't happen")
                return
            release_round_idx = r-get_num_tranches()
            seg_length=len(flexifix_constr_names)/NUM_RELEASE_ROUNDS
            
            start_idx=math.ceil(release_round_idx*seg_length)
            end_idx=math.ceil((release_round_idx+1)*seg_length)
            if release_round_idx>=NUM_RELEASE_ROUNDS-1:
                end_idx=len(flexifix_constr_names)
            
            for flexifix_name in flexifix_constr_names[start_idx:end_idx]:
                lp_problem.constraints.pop(flexifix_name)

    previously_allowed_flexi_var_names:list[str]=[]
    def set_up_flexi_variables(r:int):
        #TODO can speed up significantly by instead flagging indices of variables to enable by round (since they will always be enabled after r > some n ).
        nonlocal lp_problem
        nonlocal previously_allowed_flexi_var_names

        do_prior_solve_cut=False

        last_run_cut_var_name="lastRunCut"

        if r > 1 and do_prior_solve_cut:
            lp_problem.constraints.pop(last_run_cut_var_name)

        if KEEP_PREVIOUS_FLEXI:
            flexifix_variables(r)

        if r>=get_num_tranches():
            return

        if r > 0:
            # remove last flexi variables
            for constrs in flexible_allowed_constr_sets + flexible_forbidden_constr_sets:
                for _, constr_name in constrs: 
                    lp_problem.constraints.pop(constr_name,None) # None arg because it is okay if it does not exist
        num_allowed={k:0 for k in flexible_forbidden_constr_types}
        num_forbidden={k:0 for k in flexible_forbidden_constr_types}
        allowed_flexi_vars:list[LpVariable]=[]
        for i in range(len(flexible_allowed_constr_tranches)):
            if flexible_allowed_constr_tranches[i][r]:
                constrs_to_add=flexible_allowed_constr_sets[i]
                num_allowed[flexible_forbidden_constr_types[i]]+=1
                allowed_flexi_vars.append(flexible_vars[i])
            else:
                constrs_to_add=flexible_forbidden_constr_sets[i]
                num_forbidden[flexible_forbidden_constr_types[i]]+=1
            for constr in constrs_to_add:
                lp_problem+= constr
        
        # Forbid original geomections that are inactive and would be forbidden if . This ensures the original geomections aren't given special treatment.


        if r > 0 and do_prior_solve_cut:
            # Add constraint that lpsum(currently_active_variables) >= len(currently_active_variables)*(1 - sum(flexi_variables_allowed_for_first_time_this_round))
            # Essentially the cut that we already know from prior solve.
            new_allowed_vars=[var for var in allowed_flexi_vars if var.name not in previously_allowed_flexi_var_names]
            print("New allowed vars (first 10):",new_allowed_vars[:10])
            #print(len(allowed_flexi_vars),len(previously_allowed_flexi_vars))
            active_site_vars = []
            for site in site_var_dict:
                for from_altloc in site_var_dict[site]:
                    active_site_vars.extend(var for var in site_var_dict[site][from_altloc].values() if var.value()>0.5)
            for site in site_altposvar_dict:
                for from_altloc in site_altposvar_dict[site]:
                    active_site_vars.extend(var for var in site_altposvar_dict[site][from_altloc].values() if var.value()>0.5) 
            lp_problem+=(
                lpSum(active_site_vars)>=len(active_site_vars)*(1-lpSum(new_allowed_vars)),
                last_run_cut_var_name
            )
            
        previously_allowed_flexi_var_names=[var.name for var in allowed_flexi_vars]

        num_inactive_og_forbidden={k:0 for k in flexible_bad_original_constr_types}
        num_inactive_og_allowed={k:0 for k in flexible_bad_original_constr_types}

        for constrs in flexible_bad_original_allowed_constrs + flexible_bad_original_forbid_constrs:
            for _, constr_name in constrs: 
                lp_problem.constraints.pop(constr_name,None) # None arg because it is okay if it does not exist
        

        
        ''' # BUGGED TODO
        global pooled_method  
        def pooled_method(i):
            constrs_to_add=flexible_bad_original_allowed_constrs[i]
            tracking_tuple=None
            # Check whether the original geomection is inactive. Otherwise, do not forbid.
            if value(flexible_bad_original_vars[i]) < 0.5:
                if not always_allow_original_tranches[i][r]:
                    assert r > 0
                    tracking_tuple=(False,flexible_bad_original_constr_types[i])
                    constrs_to_add = flexible_bad_original_forbid_constrs[i]
                    #num_inactive_og_forbidden[flexible_bad_original_constr_types[i]]+=1
                else:
                    tracking_tuple=(True,flexible_bad_original_constr_types[i])
                    #num_inactive_og_allowed[flexible_bad_original_constr_types[i]]+=1
            return (constrs_to_add,tracking_tuple)
            
        with Pool(UntangleFunctions.NUM_THREADS) as p:
            list_of_constrs_to_add:list[list,type] = p.map(pooled_method,range(len(always_allow_original_tranches)))
        for constrs_to_add,tracking_tuple in list_of_constrs_to_add:
            if tracking_tuple is not None:
                if tracking_tuple[0]:
                    num_inactive_og_allowed[tracking_tuple[1]]+=1
                else:
                    num_inactive_og_forbidden[tracking_tuple[1]]+=1
            for constr in constrs_to_add:
                #assert value(constr[0]) > 0.5
                try:
                    #assert constr[1] not in lp_problem.constraints
                    lp_problem+= constr
                except Exception as e:
                    print(constr)
                    print(constrs_to_add)
                    for constrs in flexible_bad_original_allowed_constrs + flexible_bad_original_forbid_constrs:
                        for _, constr_name in constrs: 
                            if constr_name == constr[1]:
                                print("Found match (expect 1)")
                    for _,other_name in constrs_to_add:
                        if constr[1]==other_name:
                            print("Found other match (expect 1)")
                    raise e
        '''

        for i in range(len(always_allow_original_tranches)):
            constrs_to_add=flexible_bad_original_allowed_constrs[i]
            # Check whether the original geomection is inactive. Otherwise, do not forbid.
            if value(flexible_bad_original_vars[i]) < 0.5:
                if not always_allow_original_tranches[i][r]:
                    assert r > 0
                    constrs_to_add=flexible_bad_original_forbid_constrs[i]
                    num_inactive_og_forbidden[flexible_bad_original_constr_types[i]]+=1
                else:
                    num_inactive_og_allowed[flexible_bad_original_constr_types[i]]+=1

                
            for constr in constrs_to_add:
                #assert value(constr[0]) > 0.5
                try:
                    assert constr[1] not in lp_problem.constraints
                    lp_problem+= constr
                except Exception as e:
                    print(constr)
                    print(constrs_to_add)
                    for constrs in flexible_bad_original_allowed_constrs + flexible_bad_original_forbid_constrs:
                        for _, constr_name in constrs: 
                            if constr_name == constr[1]:
                                print("Found match (expect 1)")
                    for _,other_name in constrs_to_add:
                        if constr[1]==other_name:
                            print("Found other match (expect 1)")
                    raise e
                
                
         

        print(f"Num 'flexible' geomections allowed for round {r+1}: {num_allowed}")
        print(f"Num 'flexible' geomections forbidden for round {r+1}: {num_forbidden}")
        print(len(always_allow_original_tranches))
        if r > 0:
            print(f"Num inactive bad original geomections forbidden for round {r+1}: {num_inactive_og_forbidden}")
            print(f"Num inactive bad original geomections allowed for round {r+1}: {num_inactive_og_allowed}")


    def construct_conformations_from_solution():
        conf_geomection_dict:dict[str,list[LP_Input.Geomection]]={a:[] for a in all_altlocs} # Active geomections of each conformation. Keys are the conformation labels (altloc ids).
        conf_atom_chunk_dict:dict[str,list[AtomChunk]]={a:[] for a in all_altlocs}
        # smallest_unit_geomection_types = (ConstraintsHandler.BondConstraint,ConstraintsHandler.NonbondConstraint,ConstraintsHandler.ClashConstraint)

        active_geomections:list[LP_Input.Geomection]=[]
        for disordered_geomection_var_dict in mega_geomection_var_dict.values():
            active_geomections.extend([constraint for constraint,var in disordered_geomection_var_dict.values() if var.value()>0.5])
        #     for geomection in active_geomections:
        #         for other_geomection in active_geomections:
        
        

        def get_connected(chunk:AtomChunk,connected_chunks:list[AtomChunk]=[],connected_geomections:list[LP_Input.Geomection]=[]):

            connected_chunks.append(chunk)
            for geomection in active_geomections:
                if geomection in connected_geomections:
                    continue
                if chunk in geomection.atom_chunks:
                    connected_geomections.append(geomection) 
                    for other_chunk in geomection.atom_chunks:
                        if other_chunk not in connected_chunks:
                            connected_chunks,connected_geomections = get_connected(ch,connected_chunks,connected_geomections)


            return connected_chunks,connected_geomections

        '''
        for ch in chunk_sites:
            if ch.get_site_num() == lowest_site_num:
                conf_atom_chunk_dict[ch.altloc],conf_geomection_dict[ch.altloc] = get_connected(ch,active_geomections)
        '''

        for ch in chunk_sites:
            for assigned_chunks in conf_atom_chunk_dict.values():
                if ch in assigned_chunks:
                    break
            else:
                # Not assigned to any, so free to choose altloc.
                connected_chunks, connected_geomections = get_connected(ch,active_geomections)
                conf_atom_chunk_dict[ch.altloc].extend(connected_chunks)
                connected_geomections.extend(conf_geomection_dict[ch.altloc])

        return conf_atom_chunk_dict,conf_geomection_dict

    sites_restricted_by_altloc=[]
    altloc_pool=[]
    def set_up_altloc_subset_restrictions(altloc_subset_size):
        nonlocal sites_restricted_by_altloc
        nonlocal lp_problem
        nonlocal altloc_pool
        
        altloc_subset_size= min(altloc_subset_size,len(all_altlocs))

        # Free restrictions on variables from last function call
        for var_name in sites_restricted_by_altloc:
            lp_problem.constraints.pop(var_name)        
        sites_restricted_by_altloc=[]

        if altloc_subset_size==len(all_altlocs):
            return all_altlocs
        
        altlocs_to_restrict=[]
        #DEBUG_ALWAYS_HAVE_ALTLOCS=["B"] # TODO always have worst conformation. (sum up distance terms for each to_altloc)
        DEBUG_ALWAYS_HAVE_ALTLOCS=[] 
        while len(altlocs_to_restrict)<len(all_altlocs)-altloc_subset_size:
            if len(altloc_pool)==0:
                altloc_pool = [a for a in all_altlocs if a not in altlocs_to_restrict]
            if DEBUG_ALWAYS_HAVE_ALTLOCS is not None:
                altloc_pool = [a for a in altloc_pool if a not in DEBUG_ALWAYS_HAVE_ALTLOCS]
            altloc_selected = random.choice(altloc_pool)
            altloc_pool.remove(altloc_selected)
            altlocs_to_restrict.append(altloc_selected)
        for altloc, active_geomections in construct_conformations_from_solution()[1].items():
            

            for restricted_to_altloc in altlocs_to_restrict:
                var_name = get_force_no_flips_name(site,restricted_to_altloc) 
                # Note this variable name might correspond to different variables in different rounds.
                if var_name not in lp_problem.constraints:
                    for from_altloc in site_var_dict[site]:
                        if restricted_to_altloc not in site_var_dict[site][from_altloc]: 
                            continue
                        if site_var_dict[site][from_altloc][restricted_to_altloc].value()>0.5:
                            active_var= site_var_dict[site][from_altloc][restricted_to_altloc]
                            break
                    else: # This conformer site can't possibly have been assigned the restricted altloc.
                        continue
                    lp_problem += (
                        active_var==1,
                        get_force_no_flips_name(site,restricted_to_altloc)
                    ) 
                    sites_restricted_by_altloc.append(var_name)
        return [alt for alt in all_altlocs if alt not in altlocs_to_restrict]


    assert len(flexible_allowed_constr_tranches)==len(flexible_allowed_constr_sets)==len(flexible_forbidden_constr_sets), (len(flexible_allowed_constr_tranches),len(flexible_allowed_constr_sets),len(flexible_forbidden_constr_sets))
    assert len(always_allow_original_tranches)==len(flexible_bad_original_vars)==len(flexible_bad_original_allowed_constrs)==len(flexible_bad_original_forbid_constrs),\
        (len(always_allow_original_tranches),len(flexible_bad_original_vars),len(flexible_bad_original_allowed_constrs),len(flexible_bad_original_forbid_constrs))
    
    ########### LP Solver options ###########
    class Solver(Enum):
        COIN=PULP_CBC_CMD
        CPLX_PY=CPLEX_PY
        CPLX_CMD=CPLEX_CMD # NOTE Need to follow "Additional environment variables per solver" section in this guide: https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html
    pulp_solver=Solver.CPLX_CMD
    #pulp_solver = Solver.CPLX_PY
    #pulp_solver = Solver.COIN
    
    easy_timeLimit=min(DIFFICULT_SOLVE_TIME_THRESHOLD_IN_MINS,max_mins_start)
    if max_mins_start is not None:
        difficult_timeLimit=max_mins_start
    else:
        difficult_timeLimit=None
    #timeLimit=None
    logPath=log_out_dir+"solver_log.txt"
    #logPath=None

    warmStart=True
    #gapRel=0.0003
    #gapRel=0.001
    solver_class=PULP_CBC_CMD
    solver_options=[]
    if pulp_solver == Solver.COIN: 
        solver_class = PULP_CBC_CMD
    elif pulp_solver == Solver.CPLX_CMD:
        solver_class = CPLEX_CMD
    # CPLEX parameters: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/tutorials/InteractiveOptimizer/settingParams.html
    # CPLEX status: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.10.0/ilog.odms.cplex.help/refcallablelibrary/macros/Solution_status_codes.html
        #solver_options.append("set parallel -1")
        #path='~/ibm/ILOG/CPLEX_STUDIO2211/cplex'
    #https://coin-or.github.io/pulp/technical/solvers.html#pulp.apis.CPLEX_PY
    elif pulp_solver == Solver.CPLX_PY:
        solver_class = CPLEX_PY
    else:
        raise Exception("not implemented")
    #solver = solver_class(timeLimit=timeLimit,threads=THREADS,logPath=logPath,warmStart=warmStart,gapRel=gapRel,path=path)
    extra_args={}
    if pulp_solver==Solver.CPLX_CMD:
        disable_probe=False
        aggressive_probe=False
        #### Options to speed up root node relaxation: https://www.ibm.com/docs/en/icos/22.1.2?topic=problems-too-much-time-node-0
        pseudoreduced_branching=True # Y |  Computationally cheap. LEAVE THIS TRUE
        emphasize_feasibility=False  # ?
        # Other options related to root node relaxation
        detailed_display=True
        barrier_root_solve=False
        #### Options to speed up root node processing:
        turn_off_cuts=False # Y |
        disable_heuristics=True # Y |

        testing_options=True

        #####
        emphasize_optimize=False
        strong_branching=False # Computationally intensive
        diverse_solutions=False
        network_simplex=False # Employ network simplex for initial relaxation (of what? When?)
        priority_order=None # 1: decreasing cost coefficients, 2: increasing bound range, 3: increasing with matrix coefficient count 
        opportunistic_parallelism=True # Faster but not deterministic
        repeat_presolve=True
        symmetry_breaking_aggressiveness=None
        aggressive_GUB_cuts=False
        branch_up=True # We use many set partitioning constraints, which suggests this option may be good, according to Klotz & Newman 2013 DOI: 10.1016/j.sorms.2012.12.001
        max_cutting_planes_at_root=1
        disable_coefficient_reduction=False # TODO Try this
        disable_presolve=False # TODO Try this
        preprocessing_on_relaxation=True #TODO Try disabling
        #extra_args["path"]='/home/speno/ibm/ILOG/CPLEX_STUDIO2211/cplex/bin/x86-64_linux/cplex'
        extra_args["maxMemory"]=MEMORYLIMITGB*1e3
        extra_args["options"]=solver_options

        #### Perturbations https://www.ibm.com/docs/en/icos/22.1.1?topic=problems-numeric-difficulties
        always_perturb=True
        ####
        # TODO try probe
        #solver_options.append("set mip strategy probe -1") # -1: no, 0: auto, 1-3: increasingly aggressive probing
        # if len(all_altlocs) >6:
        #     solver_options.append("set mip strategy probe 3")
        # else:
        #     solver_options.append("set mip strategy probe 3")
            #solver_options.append("set mip strategy probe -1")
        # Note: Dual simplex seems to be unusually horrible for this problem. Avoid at all costs. 
        if testing_options:
            solver_options.append("set read scale 1") # aggressive scaling
            PRIMAL_SIMPLEX=1 # NOTE: Seems better than dual simplex
            BARRIER=4 # Also better than dual simplex.  # https://www.ibm.com/docs/en/icos/22.1.0?topic=performance-detecting-eliminating-dense-columns
            SIFTING=5 
            CONCURRENT=6
            
            #solver_options.append(f"set mip strategy subalgorithm {BARRIER}") # https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-mip-subproblem-algorithm
            solver_options.append(f"set mip strategy startalgorithm {PRIMAL_SIMPLEX}") # https://www.ibm.com/docs/en/icos/22.1.2?topic=problems-unsatisfactory-optimization-subproblems
            solver_options.append(f"set lpmethod {PRIMAL_SIMPLEX}") # https://www.ibm.com/docs/en/cofz/12.9.0?topic=parameters-algorithm-continuous-linear-problems
            #solver_options.append(f"set barrier limits objrange ??") # https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-barrier-objective-range
            #solver_options.append(f"set barrier algorithm 3") # https://www.ibm.com/docs/en/icos/22.1.0?topic=parameters-barrier-algorithm
            #solver_options.append("set simplex tolerances feasibility 1e-1")
            #solver_options.append("set simplex tolerances optimality 1e-1")
            #solver_options.append("set simplex pgradient -1") # Reduced cost pricing # https://www.ibm.com/docs/en/icos/22.1.1?topic=performance-simplex-parameters
            #solver_options.append("set simplex pgradient 4") # BAD
            #solver_options.append("set preprocessing dependency 3") # https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-dependency-switch
            #symmetry_breaking_aggressiveness=5
            emphasize_feasibility=True
            #aggressive_GUB_cuts=True
            #aggressive_probe=True
            max_cutting_planes_at_root=1
            turn_off_cuts=True
            disable_heuristics=True
            solver_options.append(f"set mip cuts clique 2")
            #solver_options.append(f"set mip cuts zerohalfcut 0") # Note the advice on these cuts: https://www.ibm.com/docs/en/cofz/12.10.0?topic=parameters-mip-zero-half-cuts-switch
            solver_options.append(f"set mip cuts gubcovers 2")
            #solver_options.append(f"set mip cuts gomory 0")
            #solver_options.append(f"set mip cuts implied 0") # global
            solver_options.append(f"set mip cuts localimplied 0")

            #####
            #solver_options.append(f"set mip cuts pathcut 1") # flow cover cut
            #solver_options.append(f"set mip cuts disjunctive 1")
            #solver_options.append(f"set mip cuts mircut 1")
            # TO TEST
            #branch_up=True

            solver_options.append(f"set preprocessing dual -1") # Don't solve dual https://www.ibm.com/docs/en/cofz/12.10.0?topic=parameters-presolve-dual-setting

        else: # Best-performing set of variables tested so far, which are modified above.
            solver_options.append("set read scale 1") # aggressive scaling
            PRIMAL_SIMPLEX=1 # NOTE: Seems better than dual simplex
            BARRIER=4 # Also better than dual simplex.  # https://www.ibm.com/docs/en/icos/22.1.0?topic=performance-detecting-eliminating-dense-columns
            SIFTING=5 
            CONCURRENT=6
            # Do NOT use Concurrent (or barrier?) for ???, as it will employ dual simplex after root node relaxation which can be EXTREMELY slow.
            solver_options.append(f"set mip strategy startalgorithm {PRIMAL_SIMPLEX}") # https://www.ibm.com/docs/en/icos/22.1.2?topic=problems-unsatisfactory-optimization-subproblems
            solver_options.append(f"set lpmethod {PRIMAL_SIMPLEX}") # https://www.ibm.com/docs/en/cofz/12.9.0?topic=parameters-algorithm-continuous-linear-problems
            #solver_options.append(f"set barrier limits objrange ??") # https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-barrier-objective-range
            #solver_options.append(f"set barrier algorithm 3") # https://www.ibm.com/docs/en/icos/22.1.0?topic=parameters-barrier-algorithm
            #solver_options.append("set simplex tolerances feasibility 1e-1")
            #solver_options.append("set simplex tolerances optimality 1e-1")
            #solver_options.append("set simplex pgradient -1") # Reduced cost pricing # https://www.ibm.com/docs/en/icos/22.1.1?topic=performance-simplex-parameters
            #solver_options.append("set simplex pgradient 4") # BAD
            #solver_options.append("set preprocessing dependency 3") # https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-dependency-switch
            #symmetry_breaking_aggressiveness=5
            branch_up=False
            #aggressive_GUB_cuts=True
            #aggressive_probe=True
            turn_off_cuts=True
            #max_cutting_planes_at_root=1
            emphasize_feasibility=True
            solver_options.append(f"set preprocessing dual -1") # Don't solve dual https://www.ibm.com/docs/en/cofz/12.10.0?topic=parameters-presolve-dual-setting
        if turn_off_cuts:
            solver_options.insert(0,f" set mip cuts all -1")


        if emphasize_optimize:
            solver_options.append("set emphasis mip 3") #Emphasis on moving best bound value.
        if emphasize_feasibility:
            assert not emphasize_optimize
            solver_options.append("set emphasis mip 1") 
        if disable_probe:
            solver_options.append("set mip strategy probe -1")
        if aggressive_probe:
            assert not disable_probe
            solver_options.append("set mip strategy probe 3")
            probe_time_limit_mins=5
            solver_options.append(f"set mip limits probetime {int(probe_time_limit_mins*60)}")
        #solver_options.append("set mip strategy variableselect -1")
        if strong_branching:
            solver_options.append("set mip strategy variableselect 3")
        if pseudoreduced_branching:
            assert not strong_branching
            solver_options.append("set mip strategy variableselect 4")
        if priority_order is not None:
            assert -1 <= priority_order <= 3
            solver_options.append(f"set mip ordertype {priority_order}")

        if network_simplex:
            solver_options.append(f"set mip strategy startalgorithm 3")

        if opportunistic_parallelism:
            solver_options.append("set parallel -1")


        if diverse_solutions:
            assert False, "Collecting multiple solutions this way is not yet supported"
            solver_options.append("set mip limits populate 20")
            solver_options.append("set mip pool replace 2")
        else:
            solver_options.append("set mip limits populate 1")

        if repeat_presolve and not disable_presolve:
            solver_options.append("set preprocessing repeatpresolve 3") # TODO 2 or 3?

        if symmetry_breaking_aggressiveness is not None:
            assert -1<=symmetry_breaking_aggressiveness<=5
            solver_options.append(f"set preprocessing symmetry {symmetry_breaking_aggressiveness}")
            
        if aggressive_GUB_cuts:
            solver_options.append(f"set mip cuts gubcovers 2")
        if max_cutting_planes_at_root is not None:
            solver_options.append(f"set mip limits cutpasses {max_cutting_planes_at_root}")
        if disable_heuristics:
            solver_options.append(f"set mip strategy heuristiceffort 0")
            solver_options.append(f"set mip strategy heuristicfreq -1")

        if branch_up:
            solver_options.append(f"set mip strategy branch 1")

         # https://www.ibm.com/docs/en/icos/22.1.1?topic=mip-preprocessing-presolver-aggregator
        if disable_coefficient_reduction:
            solver_options.append(f"set preprocessing coeffreduce 0")
        if disable_presolve:
            solver_options.append(f"set preprocessing presolve no")
        if detailed_display:
            solver_options.append(f"set mip display 4") # Set to 5 to show LP subproblem at nodes, and not just root.
            #solver_options.append(f"set mip display 5") 

        if (preprocessing_on_relaxation is not None) and (not disable_presolve):
            solver_options.append(f"set preprocessing relax {int(preprocessing_on_relaxation)}") # Whether to apply presolve (again/during?), to the root relaxation


        if always_perturb:
            solver_options.append(f"set simplex perturbationlimit yes 1e-06")

        assert pseudoreduced_branching
        assert not disable_probe
        #assert turn_off_cuts # TODO test which cuts are worth keeping.
        
        if barrier_root_solve:
            solver_options.append("set lpmethod 4")




        # TODO CRITICAL start from advanced basis 
        # https://www.ibm.com/docs/en/cofz/12.9.0?topic=performance-starting-from-advanced-basis
        # read bas
        # Need to modify
        # python3.12/site-packages/pulp/apis/cplex_api.py:
        # tmpLp, tmpSol, tmpMst,***tmpBas*** = self.create_tmp_files(lp.name, "lp", "sol", "mst",***"bas"***)
        # ...
        # cplex_cmds += "optimize\n"
        # cplex_cmds += "write " + tmpSol + "\n"
        # ***cplex_cmds += "write " + tmpBas + "\n" + "y\n"***
        #advanced_basis_solver_options = [s for s in solver_options]

        # advanced_basis=True
        # if advanced_basis:
        #     extra_args["keepFiles"]=True
        #     advanced_basis_solver_options.append(f"read Untangling_Problem-pulp.bas")

        #solver_options.append("set mip strategy probe -1")
    #########################################

    #print(len(flexible_allowed_constr_tranches),len(always_allow_original_tranches))
    start_time=time()


    def get_status(verbose=False):
        print("Status:", LpStatus[lp_problem.status])


        if verbose:
            for v in lp_problem.variables():
                if v.value() > 0.5:
                    print(v.name, "=", v.value())

        print(f"Target: {out_handle}")
        total_distance = value(lp_problem.objective)
        diff=total_distance/initial_badness-1
        print(f"Total distance = {total_distance} ({100*(diff):.3f}%)")
        log(f"{100*(diff):.3f}% ({total_distance})")
        #plt.scatter()


    # def good_constraint(constraint:LP_Input.Geomection):
    #     good_z_threshold=2
    #     return constraint.z_score<=good_z_threshold
    

    def get_active_to_altloc(atom_chunk:AtomChunk):
        raise NotImplementedError()
        for to_altloc, site_var in site_var_dict[VariableID.Atom(atom_chunk)][atom_chunk.altloc].items():
            if site_var.value()>0.5:
                return to_altloc
        raise Exception(f"No active to_altloc variable found for {atom_chunk.unique_id()}")
            

        
    # TODO Function that allows for getting the angle of which an atom is the centre atom of. And of which it is the endpoint of.
    # TODO Function that checks whether a bond is not possible given current .

    cement_constraints=[]
    def generate_cemented_triplet_constraints(atom_name_considered,z_bounds:tuple[float],resname_considered=None,n_to_1=True,n_to_m=True,chance_residue_considered=0.5,ignored_altlocs=[]):
        # Fix angle and its bonds to one variable if they are all good.
        nonlocal cement_constraints
        midpoint_atom_constraints=[]

        def ignore_constraint(constraint:tuple[LP_Input.Geomection,LpVariable]):
            for atom_chunk in constraint[0].atom_chunks:
                if get_active_to_altloc(atom_chunk) in ignored_altlocs:
                    return True
            return False

        for disordered_geomection_var_dict in mega_geomection_var_dict.values():
            midp_constrs = [(constraint,var) for constraint,var in disordered_geomection_var_dict.values()\
                if atom_name_considered in constraint.atom_names]
            
            midp_constrs = [midpc for midpc in midp_constrs if not ignore_constraint(midpc)]
            
            midpoint_atom_constraints.extend(midp_constrs) 
        
        good_bond_geomection_resnum_dict:dict[int,list[tuple[LP_Input.Geomection,LpVariable]]]={}
        good_angle_geomection_resnum_dict:dict[int,list[tuple[LP_Input.Geomection,LpVariable]]]={}
        if len(midpoint_atom_constraints)==0:
            return
        for constraint,var in midpoint_atom_constraints:
            if constraint.connection_type == ConstraintsHandler.BondConstraint:
                good_geomection_dict=good_bond_geomection_resnum_dict
            elif constraint.connection_type == ConstraintsHandler.AngleConstraint:
                good_geomection_dict=good_angle_geomection_resnum_dict
            else:
                continue
            if z_bounds[0]<=constraint.z_score<z_bounds[1]:
                N_idx=constraint.atom_names.index(atom_name_considered)
                resnum=constraint.res_nums[N_idx]
                mid_atm_from_altloc = constraint.from_altlocs[N_idx]
                if resnum not in good_geomection_dict:
                    good_geomection_dict[resnum]=[]
                good_geomection_dict[resnum].append((mid_atm_from_altloc,constraint,var))
        if len(good_angle_geomection_resnum_dict)==0:
            return
        #good_groups=[]
        all_cement_bond_vars:dict[str,dict[str,list[list[tuple[LP_Input.Geomection,LpVariable]]]]]={}
        for resnum, angle_constr_var_tuples in good_angle_geomection_resnum_dict.items():
            if random.random()<=chance_residue_considered:
                continue
            if resnum not in good_bond_geomection_resnum_dict:
                continue
            for (mid_frm_altl, angle_constr,angle_var) in angle_constr_var_tuples:
                atom_chunks = angle_constr.atom_chunks
                pos_option_indices=angle_constr.position_option_indices
                matching_bonds_vars=[]
                for i in range(2):
                    bond_chunks = atom_chunks[i:i+2]
                    bond_pos_indices = pos_option_indices[i:i+2]
        
                    matching=[(constr,var) for _,constr,var  in good_bond_geomection_resnum_dict[resnum] if (constr.position_option_indices==bond_pos_indices and constr.atom_chunks==bond_chunks)]
                    if len(matching)==0:
                        matching=[var for _,constr,var in good_bond_geomection_resnum_dict[resnum]
                                if (reversed(constr.position_option_indices)==bond_pos_indices and reversed(constr.atom_chunks)==bond_chunks)]
                    if len(matching)==0:
                        break
                    else:
                        assert len(matching)==1
                        matching_bonds_vars.append(matching[0])
                else:
                    #good_groups.append(matching_bonds_vars+[angle_var,])
                    #good_groups.append(matching_bonds_vars+[angle_var])
                    if resnum not in all_cement_bond_vars:
                        all_cement_bond_vars[resnum]={}
                    if mid_frm_altl not in all_cement_bond_vars[resnum]:
                        all_cement_bond_vars[resnum][mid_frm_altl]=[[] for _ in matching_bonds_vars]
                    LH_bond,RH_bond=matching_bonds_vars # !!!TODO FIXME correct order not guaranteed!!!
                    bond_lists=all_cement_bond_vars[resnum][mid_frm_altl]
                    if LH_bond not in bond_lists[0]: 
                        bond_lists[0].append(LH_bond)
                    if RH_bond not in bond_lists[1]:
                        bond_lists[1].append(RH_bond)

                    # cement_constr1=(angle_var==matching_bonds_vars[1],"cement1"+angle_var.name)
                    # cement_constr2=(matching_bonds_vars[1]==matching_bonds_vars[0],"cement2"+angle_var.name)
                    # cement_constraints.append(cement_constr1)
                    # cement_constraints.append(cement_constr2)
                    #cement_constr=(angle_var==matching_bonds_vars[1]==matching_bonds_vars[0],"cement@"+angle_var.name)
                    #cement_constraints.append(cement_constr)
        new_cement=[]
        
        ONExONE=1
        NxONE=2
        NxM=3
        cement_kind_count = {ONExONE:0,NxONE:0,NxM:0}
        for resnum in all_cement_bond_vars:
            for altloc, (bond_sets) in all_cement_bond_vars[resnum].items():
                
                iterator:itertools.combinations[tuple[list[tuple[LP_Input.Geomection,LpVariable]],list[tuple[LP_Input.Geomection,LpVariable]]]]=itertools.combinations(bond_sets, 2)
                for LH_bonds,RH_bonds in iterator:
                    # If not 1 to 1, then the rule is: if a good var on one side is active, must also choose a good var on the other side.
                    num_currently_active_LH=len([b for b in LH_bonds if b[1].value()> 0.5])
                    num_currently_active_RH=len([b for b in RH_bonds if b[1].value()>0.5])
                    #XXX
                    nameA=[name for name in LH_bonds[0][0].atom_names if (name!=atom_name_considered or len(set(LH_bonds[0][0].atom_names[0][0]))==1)][0]
                    nameB=[name for name in RH_bonds[0][0].atom_names if (name!=atom_name_considered or len(set(RH_bonds[0][0].atom_names[0][0]))==1)][0]
                    if num_currently_active_LH!=num_currently_active_RH:
                        continue # To avoid infeasibility. TODO fix by employing function that checks if forbidden. 
                                 # Then check that the max number of the sum on both sides after removing forbidden is the same
                    if len(LH_bonds) == len(RH_bonds)==1:
                        kind=ONExONE
                    elif len(LH_bonds) == 1 or len(RH_bonds)==1:
                        kind=NxONE
                    else:
                        kind=NxM
                    if not n_to_m and kind==NxM:
                        continue
                    if not n_to_1 and kind==NxONE:
                        continue
                    cement_constr=(lpSum([l[1] for l in LH_bonds])==lpSum([r[1] for r in RH_bonds]),f"cement{z_bounds[0]}to{z_bounds[1]}~{resnum}.{atom_name_considered}.{altloc}{CEMSYMBOL}{len(LH_bonds)};{len(RH_bonds)}~{nameA}@{nameB}")
                    new_cement.append(cement_constr)
                    cement_kind_count[kind]+=1

        cement_constraints.extend(new_cement)
        resname_str=f" {resname_considered}" if resname_considered is not None else ""
        for key, count in cement_kind_count.items():
            if count==0:
                continue
            print(f"Generated {cement_kind_count[key]} Z = [{z_bounds[0]}, {z_bounds[1]}) cement constraint{'s' if cement_kind_count[key]>1 else ''} of {'1:1' if key==ONExONE else 'n:1' if key==NxONE else 'n:m'} geomections for {atom_name_considered}{resname_str}")
        
    def apply_cement():
        nonlocal lp_problem
        for constr in cement_constraints:
            lp_problem+=constr
        print(f"Cemented {len(cement_constraints)} groups of geomections")
    def remove_cement():
        nonlocal lp_problem
        for constr in cement_constraints:
            lp_problem.constraints.pop(constr[1],None)

             
            
    def clear_cement_constraints():
        nonlocal cement_constraints
        cement_constraints=[]
    def generate_cement_constraints(ignored_altlocs=[]):
        # generate_cemented_triplet_constraints("C")
        # generate_cemented_triplet_constraints("N")
        # generate_cemented_triplet_constraints("CB")
        for atom_name in set(ch.name for ch in chunk_sites):
            for z_bounds in [(0,2)]:
            #for z_bounds in [(0,2),(2,4),(4,6),(6,8)]:
            # #for z_bounds in [(0,2),(1,3),(2,4),(3,5),(4,6),(5,7),(6,8)]:z
                generate_cemented_triplet_constraints(atom_name,z_bounds,ignored_altlocs=ignored_altlocs)

    peppered_fixed_geo_constr_names={}
    def pepper_fixed_geomections(freq=0.2,Z_min=0,Z_max=2,geomection_type=ConstraintsHandler.BondConstraint):
        if freq<=0:
            return
        nonlocal peppered_fixed_geo_constr_names
        nonlocal lp_problem
        
        if geomection_type not in peppered_fixed_geo_constr_names:
            peppered_fixed_geo_constr_names[geomection_type]=[]
        for disordered_geomection_var_dict in mega_geomection_var_dict.values():
            for constr,var in disordered_geomection_var_dict.values():
                if geomection_type!=constr.connection_type:
                    break
                if var.value()>0.5 and (Z_min <= constr.z_score <Z_max):
                    if random.random()<=freq:
                        constr_name="RANDOMFIX"+var.name
                        lp_problem += (
                            var==1,
                            constr_name
                        ) 
                        peppered_fixed_geo_constr_names[geomection_type].append(constr_name)
        print(f"Fixed {len(peppered_fixed_geo_constr_names[geomection_type])} {geomection_type}")

    def remove_random_fixed():
        nonlocal peppered_fixed_site_constr_names
        nonlocal peppered_fixed_geo_constr_names
        nonlocal lp_problem
        for constr_name in peppered_fixed_site_constr_names:
            lp_problem.constraints.pop(constr_name)
        peppered_fixed_site_constr_names=[]
        for geomection_type in peppered_fixed_geo_constr_names: 
            for constr_name in peppered_fixed_geo_constr_names[geomection_type]:
                lp_problem.constraints.pop(constr_name)
            peppered_fixed_geo_constr_names[geomection_type]=[]

                        
    def score_diagnostics(loop_idx,altloc_subset):
        # TODO CRITICAL split into generating data for each individual conformation. Store the information.
        return
        print("***********")
        print("Running score diagnostics")
        if reference_pdb_file is not None:
            swapped_model=get_swapped_file(reference_pdb_file,swaps_file,loop_idx)
            
            
            reference_pdb_file_subset=LP_Input.create_altloc_subset_model(reference_pdb_file,altloc_subset)
            swapped_model=LP_Input.create_altloc_subset_model(swapped_model,altloc_subset)
            
            for i, pdb in enumerate([reference_pdb_file_subset,swapped_model]):
                if i ==0:
                    print("Before")
                else:
                    print("After")
                print("========")

                if i!=0 or  havent_stored:
                    get_clashes(pdb)
                    UntangleFunctions.assess_geometry_wE(pdb)
                else:
                    pass

        print("***********")
        return

    def log_geometry_changes():
        if PLOTTING:
            all_sigma_costs:list[list[tuple[float]]]=[]
        ### Changed Connections ###
        # sigmas_i, costs_i, sigmas_f, costs_f
        changed_disordered_connections=[]
        bonds_replaced:list[LP_Input.Geomection]=[] # connections in original that are not present in solution
        #new_connections=[]
        for disordered_geomection_var_dict in mega_geomection_var_dict.values():
            # NOTE theses are lists of tuples
            active_constraints=[(constraint,var) for constraint,var in disordered_geomection_var_dict.values() if var.value()>0.5]
            original_constraints=[(constraint,var) for constraint,var in disordered_geomection_var_dict.values() if constraint.original()]

            if len(active_constraints)==len(original_constraints)==0:
                # For example, this could occur when a possible clash was avoided.
                continue

            connection_type = active_constraints[0][0].connection_type if len(active_constraints)>0 else original_constraints[0][0].connection_type
            if connection_type==ConstraintsHandler.BondConstraint:
                bonds_replaced.extend(og_constraint[0] \
                                    for og_constraint in original_constraints if og_constraint not in active_constraints)
            if PLOTTING:
                all_sigma_costs.append([(i[0].z_score,i[0].ts_distance,f[0].z_score, f[0].ts_distance) for (i,f) in zip(original_constraints,active_constraints)])

            no_change= all([constr.original() for constr,_ in active_constraints]) and len(active_constraints)==len(original_constraints)
            if no_change:
                continue

            if connection_type in [ConstraintsHandler.ClashConstraint,ConstraintsHandler.NonbondConstraint,ConstraintsHandler.TwoAtomPenalty]:
                original_constraints = [(constraint,var) for constraint,var in original_constraints if constraint.ts_distance!=0]
                active_constraints = [(constraint,var) for constraint,var in active_constraints if constraint.ts_distance!=0]

            unchanged = [ele for ele in original_constraints if ele in active_constraints]
            for ele in unchanged:
                original_constraints.remove(ele)
                active_constraints.remove(ele)
            
            if len(original_constraints)==0 and len(active_constraints)==0:
                continue

            original_cost = np.sum([constraint.ts_distance for constraint,_ in original_constraints])
            active_cost = np.sum([constraint.ts_distance for constraint,_ in active_constraints])

            original_constraints.sort(key=lambda x: x[0].ts_distance,reverse=True)
            active_constraints.sort(key=lambda x: x[0].ts_distance,reverse=True)


            #original_constraints.sort(key=lambda x: x[0].z_score,reverse=True)
            #active_constraints.sort(key=lambda x: x[0].z_score,reverse=True)

            item=(original_constraints,active_constraints,original_cost,active_cost)
            if changed_disordered_connections ==[]:
                changed_disordered_connections=[item]
            else:
                changed_disordered_connections.append(item)
        bonds_replaced.sort(key=lambda x: x.atom_chunks[0].get_site_num())
        bonds_replaced_each_loop.append(bonds_replaced)
        changed_disordered_connections.sort(key=lambda x: x[2]-x[3],reverse=True)              
        out_str=""
        total_distance = value(lp_problem.objective)
        diff=total_distance/initial_badness-1
        out_str+=f"Total distance = {total_distance} ({100*(diff):.3f}%)\n"
        out_str+="name, z-score, cost\n"
        for original_constraints,active_constraints,original_cost,active_cost in changed_disordered_connections:
            Dordered_constr_ref=active_constraints if len(active_constraints)>0 else original_constraints
            out_str+=Dordered_constr_ref[0][0].get_disordered_connection_id()+"\n"
            def add_disordered_block_to_str(constraints_vars:list[tuple[LP_Input.Geomection,LpVariable]]):
                nonlocal out_str
                for constraint, _var  in constraints_vars:
                    poschange_str= f"position_changes{constraint.poschange_tag}" if constraint.involves_position_changes() else ""
                    altloc_str = ','.join(constraint.from_altlocs)
                    out_str+=f"{altloc_str} {constraint.z_score:.2e} {constraint.ts_distance:.2e} {poschange_str}\n"
            
            cost_str = ""
            if original_cost>0:
                cost_str = f" {(active_cost/original_cost-1)*100:.2f}% |"
            out_str+=f"Change:{cost_str} {original_cost:.2e} --> {active_cost:.2e}\n"
            out_str+="Original\n"
            add_disordered_block_to_str(original_constraints)
            out_str+="Active\n"
            add_disordered_block_to_str(active_constraints)
            out_str+="***********************\n"
        with open(f"{log_out_dir}/ChangedConnections-{l+1}.txt",'w') as f:
            f.write(out_str)

        def write_bonds_replaced_to_file(bonds_replaced:list[LP_Input.Geomection],out_path):
            out_str="Replaced bonds\n"
            out_str+=f"Total distance = {total_distance} ({100*(diff):.3f}%)\n"
            out_str+="==============\n"
            
            bond_altloc_dict={}
            for bond in bonds_replaced:
                key = bond.get_disordered_connection_id()
                if key not in bond_altloc_dict:
                    bond_altloc_dict[key]=""
                assert bond.original()
                bond_altloc_dict[key]+=bond.from_altlocs[-1] #FIXME assumes that child altlocs will only ever be second element of from_altlocs.
            for disordered_id, altlocs in bond_altloc_dict.items():
                assert len(bond.atom_chunks)==2
                out_str+= f"{disordered_id} {altlocs}\n"
            with open(out_path,'w') as f:
                f.write(out_str)
        def write_poschanges_to_file(out_path):
            out_str="Changed coords\n"
            out_str+=f"Total distance = {total_distance} ({100*(diff):.3f}%)\n"
            out_str+="==============\n"
            site_altloc_dict={}
            for site in site_altposvar_dict:
                for from_altloc in site_altposvar_dict[site]:
                    if site_altposvar_dict[site][from_altloc][0].value()<0.5:
                        if site not in site_altloc_dict:
                            site_altloc_dict[site]=""
                        site_altloc_dict[site]+=from_altloc
            if len(site_altposvar_dict)==0:
                return # No changes
            for site, altlocs in site_altloc_dict.items():
                out_str+= f"{site} {altlocs}\n"
            with open(out_path,'w') as f:
                f.write(out_str)
        write_bonds_replaced_to_file(bonds_replaced,f"{log_out_dir}/ChangedBonds-{l+1}.txt")

        if len(site_altposvar_dict)>0:
            write_poschanges_to_file(f"{log_out_dir}/ChangedCoords-{l+1}.txt")


        if PLOTTING:
            try:
                all_sigma_costs = np.array(all_sigma_costs,dtype=np.float32)
                xlim_dict:dict[str,tuple[float]]={}
                for i, name in enumerate(["z-score_i","costs_i","z-score_f","costs_f"]):
                    X=all_sigma_costs[...,i].flatten()
                    # Same x limits for initial and final (TODO same y limits... need to refactor)
                    variable_kind = name.split("_")[0] #XXX
                    if variable_kind not in xlim_dict:
                        xlim=(np.quantile(X,0.95),np.max(X))
                        xlim_dict[variable_kind]=xlim
                    else:
                        xlim=xlim_dict[variable_kind]
                        
                    matplotlib.use("Agg")
                    plt.hist(X,bins=20,range=xlim)
                    plt.yscale('log')
                    plt.ylim([0.9,None])

                    plt.xlabel(name)
                    plt.ylabel("frequency")
                    plt.savefig(f"{name}.png")
                    plt.close()
            except Exception as e:
                print(f"Plotting failed. Error: {e}")    
    ########

    get_clashes(reference_pdb_file)
    for l in range(num_solutions):
        if l > 0 and l <= len(forced_swap_solutions):
            lp_problem.constraints.pop("forcedSwap")
        if l < len(forced_swap_solutions):
            lp_problem += (
                forced_swap_solutions[l],
                f"forcedSwap"
            )
        
        print()
        print(f"------- Start loop {l+1} of {num_solutions} ---------")
        if l == 0:
            print("Solving best solution")
        elif l == 1:
            print("Solving next-best solution")
        else:
            print(f"Solving {l+1}th-best solution")
        print(f"-----------------------------------------------------")
        print()

        log(f"Loop {l+1}/{num_solutions}")
        ## Initialise variables that will update every few rounds, for use in updating the swaps file containing the solutions. 
        site_assignment_arrays.append({})
        distances.append(np.inf)
        ##

        num_rounds=get_num_rounds()
        difficult=False
        
        for r in range(num_rounds):

            altloc_subset_sizes=[None]
            if not difficult and (ALTLOC_RUN_SUBSET_SIZES is not None):
                altloc_subset_sizes=ALTLOC_RUN_SUBSET_SIZES
            elif difficult and (ALTLOC_RUN_SUBSET_SIZES_AFTER_DIFFICULT is not None):
                altloc_subset_sizes=ALTLOC_RUN_SUBSET_SIZES_AFTER_DIFFICULT
            num_altloc_subset_runs=len(altloc_subset_sizes)

            if num_rounds>1:
                print(f"Round {r+1}/{num_rounds}")
                log(f"Round {r+1}/{num_rounds}")
                if r < len(improvement_factors_to_tolerate):
                    print(improvement_factors_to_tolerate[r])
                    log(f"{improvement_factors_to_tolerate[r]}")

            set_up_flexi_variables(r)


            for j in range(num_altloc_subset_runs):
                if altloc_subset_sizes[j] is not None and altloc_subset_sizes[j]<len(all_altlocs):
                    altlocs_in_problem = set_up_altloc_subset_restrictions(altloc_subset_sizes[j])
                    restricted_text=f"Altlocs restricted to {altlocs_in_problem}"
                    is_subset_run = False
                else:
                    restricted_text = "No altlocs restricted"
                    altlocs_in_problem=set_up_altloc_subset_restrictions(np.inf)
                    is_subset_run = True
                print(f"{restricted_text} ({j+1}/{num_altloc_subset_runs})")

                remove_cement()
                if len(altlocs_in_problem)>= MIN_ALTLOCS_TO_GLUE_GOOD_GEOMETRY_GROUPS:
                    clear_cement_constraints()
                    generate_cement_constraints(ignored_altlocs = [a for a in all_altlocs if a not in altlocs_in_problem])
                    apply_cement()

                remove_random_fixed()
                if len(altlocs_in_problem)>=MIN_ALTLOCS_TO_FIX_RANDOM:
                    pepper_fixed_sites()
                    pepper_fixed_geomections(0.2,Z_min=0,Z_max=2)
                    pepper_fixed_geomections(0.1,Z_min=2,Z_max=3)
                                
                def run_solve():
                    if create_initial_variable_files:
                        ###
                        lp_problem.writeLP(f"{log_out_dir}/LP.lp")    
                        ###
                        with open(f"{log_out_dir}/ProblemStatusStart.txt",'w') as f:
                            f.write(f"Status: {LpStatus[lp_problem.status]}\n")

                            for v in lp_problem.variables():
                                try:
                                    if v.value() > 0.5:
                                        f.write(f"{v.name} = {v.value()}\n")
                                except:
                                    raise Exception(v,v.value())
                            f.write(f"Total distance = {value(lp_problem.objective)}")
                    minutes = easy_timeLimit if not difficult else difficult_timeLimit
                    # if is_subset_run:
                    #     minutes= max(1,minutes/4)

                    # if pulp_solver==Solver.CPLX_CMD and advanced_basis:
                    #     if l==r==j==0 or os.path.getsize("Untangling_Problem-pulp.bas")<=1000:
                    #         extra_args["options"]=solver_options
                    #     else:
                    #         extra_args["options"]=advanced_basis_solver_options
                    
                    # TODO try reusing solver object and just modify member variables. Might be that pulp naturally supports loading file information from previous CPLEX run, as CPLEX environment naturally does. 
                    solver = solver_class(timeLimit=60*minutes,threads=THREADS,warmStart=warmStart,logPath=logPath,
                                          gapRel=gapRel,
                    **extra_args)
                    run_start_time=time()
                    lp_problem.solve(solver)
                    solve_time = time()-run_start_time
                    sleep(1)
                    gc.collect()

                    return solve_time
                
                solve_time = run_solve()
                if not difficult and solve_time/60>=DIFFICULT_SOLVE_TIME_THRESHOLD_IN_MINS:
                    print("shifting gears to difficult problem mode")
                    difficult=True
                    solve_time=run_solve()
                total_solve_time=time()-start_time
                log(f"Solver time: {int(solve_time/60)} m {int(solve_time%60)} s, Total: {int(total_solve_time/60)} m {int(total_solve_time%60)} s")

                print()
                print("Solver finished")
                # if solve_time/60>=max_mins*0.95:
                #     pass

                get_status(verbose=False)




                # dry = None
                # for val in connections.values():
                #     for v in val.values():
                #         dry=v.calculated_dry()
                #         break
                # tag = "_dry" if dry else ""
                with open(f"{log_out_dir}/ProblemStatusEnd.txt",'w') as f:
                    f.write(f"Status: {LpStatus[lp_problem.status]}\n")

                    for v in lp_problem.variables():
                        try:
                            if v.value() > 0.5 or v.name.startswith("Debug"):
                                f.write(f"{v.name} = {v.value()}\n")
                        except:
                            raise Exception(v,v.value())
                    f.write(f"Total distance = {value(lp_problem.objective)}")

                if lp_problem.sol_status==LpStatusInfeasible:
                    assert l>0, "Solution was infeasible!"

                    print()
                    print(f"WARNING: Finding solution {l+1} on round {r} was infeasible!")
                    break
                ##Write Active Connections##
                write_current_connections(f"{log_out_dir}/ActiveConnections.txt")
                ############################

                site_assignments:dict[VariableID,dict[str,dict[str,int]]] = {}
                site_assignment_arrays[-1]=site_assignments
                distances[-1]= value(lp_problem.objective)

                # Determine which atom has been assigned where.
                for site in site_var_dict:
                    site_assignments[site]={}
                    for from_altloc in site_var_dict[site]:
                        if lp_problem.sol_status==LpStatusInfeasible:
                            assert False
                            site_assignments[site][from_altloc]=from_altloc 
                            break 
                        poschange_str= ""
                        if site in site_altposvar_dict and from_altloc in site_altposvar_dict[site]:
                            for poschange_index, altposvar in site_altposvar_dict[site][from_altloc].items():
                                if poschange_index!=0 and altposvar.value()>0.5:
                                    assert poschange_str == ""
                                    poschange_str = f" new_position={site_altpos_dict[site][from_altloc][poschange_index].get_coord()}" 
                        to_altloc_found = False
                        for to_altloc in site_var_dict[site][from_altloc]:
                            if site_var_dict[site][from_altloc][to_altloc].value()>0.5:  # For some reason CPLEX outptuts values like 1.0000000000094025 sometimes.
                                assert not to_altloc_found
                                to_altloc_found=True
                                site_assignments[site][from_altloc]=to_altloc+poschange_str

                                            
                        assert to_altloc_found, (site, from_altloc)
                
                update_swaps_file(distances,site_assignment_arrays)  #,record_notable_improvements_threshold=0.03)
                
                log_geometry_changes()
                score_diagnostics(l,altlocs_in_problem)



            
        ### lth solution found ####
            
        if lp_problem.sol_status==LpStatusInfeasible:
            assert l>0, "Solution was infeasible!"

            print()
            print(f"WARNING: Finding solution {l+1} was infeasible! Ending solution search")
            break

       


                



        #  Conditions  imposed for next best solutions
        include_alt_pos_as_next_best_options=False
        flip_variables=[]
        if not forbid_solutions_composed_of_better_solutions:
            inverted_active_variables = []
            inverted_active_variable_names = []
        else:
            assert False


        include_alt_pos_as_next_best_options=False # Allow solutions that only differ by alternate positions # TODO Need to forbid solutions that differ only by alternate positions and by a label reassignment (it is equivalent to no label reassignment as position changes are blind to conformer labels) 
        
        REQUIRE_ONE_UNUSED_BOND_GEOMECTION=True
        
        if not REQUIRE_ONE_UNUSED_BOND_GEOMECTION:
            for chunk in chunk_sites:
                site = VariableID.Atom(chunk)
                if not site_being_considered(site):
                    continue
                if protein_sites and site.is_water: # not interested in different solutions for water. 
                    continue # This means water atoms may or may not swap for the single solution where no protein atoms swap.
                from_altloc = chunk.get_altloc()

                vars_to_force_change = list(site_var_dict[site][from_altloc].values())
                if site in site_altposvar_dict and include_alt_pos_as_next_best_options:
                    vars_to_force_change.extend(site_altposvar_dict[site][from_altloc].values())
                
                
                for var in vars_to_force_change:
                    val = var.value()
                    flip_variables.append(var)
                    if (val > 0.5) and not (str(var) in inverted_active_variable_names):
                        inverted_active_variables.append(1-var)
                        inverted_active_variable_names.append(str(var))


            if len(inverted_active_variables) == 0:
                # Require one variable to be flipped
                lp_problem += pulp.lpSum(flip_variables) >= 1, f"force_swaps_loop_{l}"  # TODO remove?
            else:
                # require at least one assignment to be different
                lp_problem += pulp.lpSum(inverted_active_variables) >= 1, f"force_next_best_solution_{l}"
            

        else: 
            geom_vars_to_force_change=[]
            active_bond_vars=[]
            for disordered_geomection_var_dict in mega_geomection_var_dict.values():
                if len(disordered_geomection_var_dict)==0:
                    continue
                reference_conn=list(disordered_geomection_var_dict.values())[0][0]
                if reference_conn.connection_type!=ConstraintsHandler.BondConstraint:
                    continue
                geom_vars_to_force_change.extend(var for _, var in disordered_geomection_var_dict.values())
                active_bond_vars.extend([var for constraint,var in disordered_geomection_var_dict.values() if var.value()>0.5])
            #active_altpos_vars=[]
            # for site in site_altposvar_dict:
            #     for from_altloc in site_altposvar_dict[site]:
            #         geom_vars_to_force_change.extend([var for var in site_altposvar_dict[site][from_altloc].values()])
            #         active_altpos_vars.extend([var for var in site_altposvar_dict[site][from_altloc].values() if var.value()>0.5])
            new_vars=[]
            #for var in (active_bond_vars+active_altpos_vars):
            for var in (active_bond_vars):
                if not (str(var) in previously_used_variable_names):
                    new_vars.append(str(var))
            previously_used_variable_names.extend(new_vars)
            if l>0:
                #print("New vars (bonds/positions):",new_vars)
                print("New bond geomections:",new_vars)
            unused_variables = [var for var in geom_vars_to_force_change if str(var) not in previously_used_variable_names]
            lp_problem += pulp.lpSum(unused_variables) >= 1, f"force_unused_bond_{l}"

        FORBID_UNUSED_O_ALTPOS=True  # TODO should forbid only if there is no change in other bond geomections with C 
        if FORBID_UNUSED_O_ALTPOS and l==0:
            active_O_altpos_vars=[]
            for site in site_altposvar_dict:
                if site.atom_name!="O":
                    continue
                for from_altloc in site_altposvar_dict[site]:
                    active_O_altpos_vars.extend([var for var in site_altposvar_dict[site][from_altloc].values() if var.value()>0.5])
            lp_problem += pulp.lpSum(active_O_altpos_vars) == len(active_O_altpos_vars), f"forbid_unused_O_altpos_{l}"


        '''
        non_original_variables=[]
        non_original_variable_history.append(non_original_variables)
        
        for disordered_geomection_var_dict in mega_geomection_var_dict.values():
            active_bonds=[(constraint,var) for constraint,var in disordered_geomection_var_dict.values() if var.value()>0.5]
            reference_conn=active_bonds[0][0]
            if reference_conn.connection_type!=ConstraintsHandler.BondConstraint:
                continue

            

            for constraint,var in active_bonds:
                val = var.value()
                flip_variables.append(var)
                if (val > 0.5): 
                    if not constraint.original():
                        non_original_variables.append(var)
                    if not (str(var) in inverted_active_variable_names):
                        inverted_active_variables.append(1-var)
                        inverted_active_variable_names.append(str(var))

            if not include_alt_pos_as_next_best_options:
                assert False
                for constraint in [constraint for constraint,var in active_bonds]:
                    if not constraint.involves_position_changes():
                        continue
                    also_exclude = [var for constraint,var in disordered_geomection_var_dict.values() if\
                                    constraint.position_option_indices==reference_conn.position_option_indices\
                                    and var.value()<0.5]
                    for var in also_exclude:
                        val = var.value()
                        flip_variables.append(var)
                        if not (str(var) in inverted_active_variable_names):
                            inverted_active_variables.append(1-var)
                            inverted_active_variable_names.append(str(var))
    
        # TODO unions of more than two ?
        for j, previous_non_original_vars in enumerate(non_original_variable_history[:-1]):
            combined_list = []
            for var in non_original_variables + previous_non_original_vars:
                if var not in combined_list:
                    combined_list.append(var)
                
                combined_list.append(all_mutually_exclusive_vars)
            lp_problem += pulp.lpSum(combined_list) >= 1, f"forbid_union_loop_{l}_idx_{j}"  # TODO remove?
        '''

        if max_bond_changes_tuple is not None:
            assert MAX_BOND_CHANGES_SECOND_HALF_ONLY
            assert l <= int(num_solutions/2) 
            if l == int(num_solutions/2):
                limit_bond_changes()





    del lp_problem
    gc.collect()
        ##################

    assert swaps_file==swaps_file_path(out_dir,out_handle,all_altlocs) # XXX
    return swaps_file,bonds_replaced_each_loop, distances


            
def get_swapped_file(unswapped_pdb_file,swap_file_path,swap_idx):

    swapper = Swapper()
    swapper.add_candidates(swap_file_path) #

    i=0
    while swapper.solutions_remaining()>0: # XXX TODO swapper should be updated to facilitate swapping the model at the specified idx.
        swapped_model, swapGroup = swapper.run(unswapped_pdb_file)
        if i==swap_idx:
            break
        i+=1
    else:
        print(f"get_clashes error - did not get model at index {swap_idx}")
        return
    return swapped_model
    
def get_clashes(pdb_model):
    clash_score_program= os.path.join(UntangleFunctions.UNTANGLER_WORKING_DIRECTORY,"Measures","clash_score_keepH.sh")

    args=["bash", clash_score_program, pdb_model]
    print (f"|+ Running: {' '.join(args)}")
    subprocess.run(args)
    

def swaps_file_path(out_dir,out_handle,altlocs):
    return f"{out_dir}/xLO-toFlip_{out_handle}-{''.join(sorted(altlocs))}.json"

if __name__=="__main__":
    handle = sys.argv[1]
    out_dir = f"{os.path.abspath(os.getcwd())}/output/{handle}/"
    os.makedirs(out_dir,exist_ok=True)
    pdbPath = f"{UntangleFunctions.pdb_data_dir()}/{handle}.pdb"
    finest_chunks,disordered_connections= LP_Input(pdbPath).calculate_paths(
        atoms_only=True,
    )
    solve(finest_chunks,disordered_connections,handle)

