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

DEBUG_FIRST_100_SITES=False


import pulp as pl
import numpy as np
from enum import Enum
from LinearOptimizer.Input import *
import itertools
import UntangleFunctions
import json
from copy import deepcopy
import gc; 
#import pulp as pl
# just for residues for now


#def solve(chunk_sites: list[Chunk],connections:dict[str,dict[str,MTSP_Solver.ChunkConnection]],out_handle:str): # 
def solve(chunk_sites: dict[str,AtomChunk],disordered_connections:dict[str,list[MTSP_Solver.AtomChunkConnection]],out_dir,out_handle:str,force_no_flips=False,num_solutions=20,force_sulfur_bridge_swap_solutions=True,
          inert_protein_sites=False,protein_sites:bool=True,water_sites:bool=True,max_mins_start=5,mins_extra_per_loop=0,gapRel=0.001): # 
    # protein_sites, water_sites: Whether these can be swapped.
    # gaprel : relative gap tolerance for the solver to stop (fraction) # https://coin-or.github.io/pulp/technical/solvers.html
    if gapRel == 0:
        gapRel = None

    print("Initialising solver")

    num_forbidden_connections=0
    num_allowed_connections=0

    if inert_protein_sites:
        assert protein_sites
    #first_atom_id = "1.N"
    lowest_site_num = np.inf
    for chunk_site in chunk_sites.values():
        if chunk_site.get_site_num() < lowest_site_num:
            lowest_site_num = chunk_site.get_site_num()


    #nodes = [chunk.unique_id() for chunk in chunk_sites.values()]

    def get_variables(unique_id)->dict:
        d={}
        d["depth"] = unique_id.split("&")[0]
        d["site_num"],d["altloc"] = unique_id.split("&")[1].split(".")
        d["name"] = unique_id.split("&")[2]
        return d
    
    def get(unique_id,key):
        return get_variables(unique_id)[key]

    def site_being_considered(site:'VariableID'):
       if DEBUG_FIRST_100_SITES and site.site_num>100:
           return False
       return (((protein_sites and not site.is_water)
             or (water_sites and site.is_water ))
             )
             #and site.site_num < 10) #XXX DEBUG TEMPORARY   


    forced_swap_solutions=[]




    lp_problem = pl.LpProblem("Untangling_Problem", LpMinimize)

    

   

    ##### Constraints ####


    class VariableKind(Enum):
        Atom="Atom"
        Bond = "Bond"
        Nonbond="Nonbond"
        Clash="Clash"
        Angle = "Angle"

    #Disordered variable id
    # Messy...
    class VariableID:
        @staticmethod
        def Atom(chunk:AtomChunk):
            return VariableID(
                chunk.get_disordered_tag(), #f"{chunk.resnum}.{chunk.name}",
                VariableKind.Atom,
                chunk.get_site_num(),
                chunk.is_water)
        def __init__(self,name:str,kind:VariableKind,site_num=None,is_water=None):
            self.name=str(name)
            self.kind = kind
            self.site_num=site_num
            self.is_water = is_water
        def __repr__(self):
            return self.name
        def __hash__(self):
            return hash((self.name, self.kind))
        def __eq__(self,other):
            return (self.name,self.kind) == (other.name,other.kind) 
        def __ne__(self,other):
            return not(self == other) 
    #bond_choices = {} # Each must sum to 1
    distance_vars = [] 
    # NOTE When assigning LP variable names, the "class" of variable should follow format of variableName_other_stuff   (class of variable is variable_type.split("_")[0]) 
    constraint_var_dict:dict[VariableID,tuple[MTSP_Solver.AtomChunkConnection,pl.LpVariable]] = {} 
    site_var_dict:dict[VariableID,dict[str,dict[str,pl.LpVariable]]] ={}
    var_dictionaries = dict(
        constraints = constraint_var_dict,
        sites = site_var_dict,
    )

    # Setup where atoms are assigned.
    # Track which are which based on original altloc.

    all_altlocs:list[str]=[]
    for ch in chunk_sites.values():
        if ch.altloc not in all_altlocs:
            all_altlocs.append(ch.altloc)
    # def site_variable_name(site: VariableID, from_altloc:str, to_altloc:str):
    #     return f"atomState_{site}_{from_altloc}.{to_altloc}"
    disordered_atom_sites:list[VariableID] = []
    
    # Setup atom swap variables
    #Crucial TODO: Too many swap variables. E.g. with two altlocs we will have four variables per atom site when we only need one.
    for i, (atom_id, chunk) in enumerate(chunk_sites.items()):
        site = VariableID.Atom(chunk)

        if not site_being_considered(site):
            continue

        if site not in disordered_atom_sites:
            disordered_atom_sites.append(site)

        from_altloc = get(chunk.unique_id(),"altloc")
        
        #TODO!!!!!!! if forbidden due to being absent from any disordered connections that involve the site, exclude! 
        if site not in site_var_dict:
            site_var_dict[site] = {}
        if from_altloc not in site_var_dict[site]:
            site_var_dict[site][from_altloc]={}


    # TODO make these variables binaries that correspond to every permutation (for more than 2 altlocs). E.g. 6 variables for 3 altlocs (could try 3)    
    dummy_one = pl.LpVariable("One",1,1,cat=const.LpBinary)
    dummy_one.setInitialValue(1)

    
    for site in disordered_atom_sites:
        if not site_being_considered(site) :
            continue
        site_altlocs = []
        for possible_altloc in site_var_dict[site]:
            site_altlocs.append(possible_altloc)


        if len(all_altlocs)>2: #TODO optimize
            for from_altloc in site_altlocs:
                # Create variable for each possible swap to other altloc
                for to_altloc in all_altlocs:
                    var_atom_assignment =  pl.LpVariable(
                        f"{site}_{from_altloc}.{to_altloc}",
                        lowBound=0,upBound=1,cat=pl.LpBinary #TODO pl.LpBinary
                    )
                    # For warm start.
                    var_atom_assignment.setInitialValue(0)
                    if to_altloc == from_altloc: 
                        var_atom_assignment.setInitialValue(1)
                    site_var_dict[site][from_altloc][to_altloc]=var_atom_assignment

                # Each ordered atom is assigned to one conformer from:to == n:1
                lp_problem += (  
                    #lpSum(site_var_dict[site][from_altloc])==1,
                    lpSum(site_var_dict[site][from_altloc].values())==1,
                    f"fromAltLoc_{site}_{from_altloc}"
                )

            # Each conformer is assigned one ordered atom. from:to == 1:n
            for to_altloc in all_altlocs:
                to_altloc_vars:list[LpVariable] = []
                for from_alt_loc_dict in site_var_dict[site].values():
                    to_altloc_vars.append(from_alt_loc_dict[to_altloc])
                lp_problem += (  
                    lpSum(to_altloc_vars)==1,
                    f"toAltLoc_{site}.{to_altloc}"
                )
        else:
            var_flipped =  pl.LpVariable(
                f"{site}_Flipped",
                lowBound=0,upBound=1,cat=pl.LpBinary #TODO pl.LpBinary
            )
            var_flipped.setInitialValue(0)
            var_not_flipped = 1-var_flipped 

            if len(site_altlocs)==1:  
                site_var_dict[site][site_altlocs[0]][site_altlocs[0]] = var_not_flipped
                other_altlocs = [a for a in all_altlocs if a != site_altlocs[0]]
                assert len(other_altlocs)==1
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
        if (force_no_flips
            or (site.site_num == lowest_site_num) # Anchor solution to one where first disordered atom is unchanged.  
            or (not site.is_water and inert_protein_sites)
        ) :
            
            for altloc in site_altlocs:
                lp_problem += (
                    site_var_dict[site][altloc][altloc]==1,
                    f"forceNoFlips_{site}_{altloc}"
                )  


    # TODO improve terminology
    # A "connection" just refers to a group of atoms with a constraint assigned by LinearOptimizer.Input.
    # A disordered connection refers to all the *possible* groupings of these atoms.

    def add_constraints_from_disordered_connection(constraint_type:VariableKind,disordered_connection: list[MTSP_Solver.AtomChunkConnection],global_score_tolerate_threshold=0):
        # Rule: If all atom assignments corresponding to a connection are active,
        # then all those atoms must be swapped to the same assignment.
        nonlocal lp_problem  



        altlocs = None
        all_combos=True

        site_altlocs_same=True
        for ch in disordered_connection[0].atom_chunks:
            site = VariableID.Atom(ch)
            if not site_being_considered(site):
                    return
            if altlocs is None:
                altlocs = set(site_var_dict[site].keys())
            else:
                if altlocs != set(site_var_dict[site].keys()):
                    site_altlocs_same=False
                    #all_combos=False
                    #assert False, "Missing altloc option for site"
                    #print(f"Warning: altlocs don't match, skipping {disordered_connection[0].get_disordered_connection_id()}")
                    #return

        # dicts indexed by code corresponding to from altlocs (e.g. "ACB" means connecting up site 1 altloc A, site 2 altloc C, site 3 altloc B)
        connection_dict:dict[str,MTSP_Solver.AtomChunkConnection]={}  
        connection_var_dict={} # indexed by from_altloc
        if site_altlocs_same:
            n=len(altlocs) 
            m=len(disordered_connection[0].atom_chunks)

            for ordered_connection_option in disordered_connection:
                altlocs_key = ''.join(ordered_connection_option.from_altlocs)
                connection_dict[altlocs_key]=ordered_connection_option

            all_combos = (len(connection_dict)==n**m)
        else:
            all_combos=False

        #assert all_combos, (disordered_connection[0].connection_type, len(connection_dict),n**m,n,m)


        for ordered_connection_option in disordered_connection:
            from_ordered_atoms = "|".join([f"{ch.resnum}.{ch.name}_{ch.altloc}" for ch in ordered_connection_option.atom_chunks])
            tag=from_ordered_atoms
            extra_tag=""
            if ordered_connection_option.hydrogen_tag!="":
                extra_tag = "Htag["+ordered_connection_option.hydrogen_tag+"]_"
            tag+=extra_tag
            ########

            # if not allowed: 
            #     for to_altloc in all_altlocs:
            #         assignment_vars = [site_var_dict[VariableID.Atom(ch)][ch.get_altloc()][to_altloc] for ch in ordered_connection_option.atom_chunks]
            #         lp_problem += (
            #             lpSum(assignment_vars) <=  len(assignment_vars)-1,   
            #             f"FORBID{constraint_type.value}_{tag}>>{to_altloc}"
            #             )
            #else:

            var_active = pl.LpVariable(f"{constraint_type.value}_{tag}",  #TODO cat=pl.LpBinary
                                lowBound=0,upBound=1,cat=pl.LpBinary)
            
            var_active.setInitialValue(0) 
            if len(set([ch.get_altloc() for ch in ordered_connection_option.atom_chunks])) == 1:
                var_active.setInitialValue(1)
            constraint_var_dict[VariableID(from_ordered_atoms+extra_tag,constraint_type.value)]=(ordered_connection_option,var_active)
            #group_vars.append(var_active)
            altlocs_key = ''.join(ordered_connection_option.from_altlocs)
            connection_var_dict[altlocs_key]=var_active



        # Get pairs corresponding to swap vs not swap.
        allowed_connections = []
        group_choices = None
        all_allowed=False
        #if (len(altlocs)==2 or debugging_more_than_2) and (disordered_connection[0].connection_type in [VariableKind.Bond.value,VariableKind.Angle.value]):
        #if not all_allowed and (len(altlocs)==2 or debugging_more_than_2) and (disordered_connection[0].connection_type!=VariableKind.Clash.value):
        
        # Find all possible groups of ordered connections that can be chosen as a choice for this disorderd connection 
        # (e.g. 3 conformers, there will be three ordered connections). Forbid really bad choices
        # TODO should still run if all_combos==False. Just forbid groups that are incomplete.
        """
        if all_combos and site_altlocs_same and not all_allowed:
            # Skip clashes since may not exist for all possible group_choices.

            #TODO test with more than 2 altlocs!
            combos:List[List[MTSP_Solver.AtomChunkConnection]] = []
            connections_assigned_to_permut=[]


            # Construct a code that corresponds to recipe for a possible choice for this connection.
            # For example, ACB_BCA would corespond to ordered connections A-B,C-C,B-A. (This would mean swapping the altlocs of one of the two A and B atoms, but leving the C altlocs unchanged)  
            # TODO Set this once per solve(), not each time constraint added
            possible_site_altlocs=[""]
            for j in range(m): # for each site
                new_possible_site_altlocs=[]
                for label in possible_site_altlocs: 
                    if label!="":
                        label+="_"

                    # Break symmetry by making group 1 correspond to first position having first altloc, group 2 having second altloc at first position, and so on.
                    if j == 0:
                        permutations = itertools.combinations(altlocs,len(altlocs)) # == [altlocs]
                    else:
                        permutations = itertools.permutations(altlocs,len(altlocs)) 
                    for permutation in permutations:
                        disordered_site=''.join(permutation) # The order the from altlocs are taken from by each group. 
                        disordered_sequence = label+disordered_site 
                        new_possible_site_altlocs.append(disordered_sequence)
                possible_site_altlocs=new_possible_site_altlocs

            group_choices:list[list[MTSP_Solver.AtomChunkConnection]]=[]
            group_choice_codes=[]
            group_choice_var_sets=[]
            for code in possible_site_altlocs:
                group=[]
                var_group=[]
                disordered_site_codes = code.split("_") # "ACB_BCA"=>["ACB", "BCA"]
                for k in range(n): # for each to_altloc (i.e. each ordered connection in this group) 
                    site_altlocs=''.join([site_codes[k] for site_codes in disordered_site_codes]) # ["ACB", "BCA"]  k=0:  => site_altlocs=AB
                    group.append(connection_dict[site_altlocs])
                    var_group.append(connection_var_dict[site_altlocs])
                group_choices.append(group)
                group_choice_var_sets.append(var_group)
                group_choice_codes.append(code)

            # Number of groups we expect
            '''
            num_groups=1
            k=n
            while k >0:
                num_groups*=k**m
                k-=1
            num_groups/=math.factorial(n) # Group order doesn't matter
            num_groups=int(num_groups)
            assert len(group_choices)==num_groups, (disordered_connection[0].connection_type, len(group_choices),num_groups, '\n'.join([str([c.from_altlocs for c in p]) for p in group_choices]),disordered_connection[0].get_disordered_connection_id())
            '''

            '''
            perm_scores=[]
            sorted_perm_scores=[] #XXX
            for group in group_choices:
                score = sum(conn.ts_distance for conn in group )
                perm_scores.append(score)
                sorted_perm_scores.append(score)
            sorted_perm_scores.sort()
            best_score,second_best_score = sorted_perm_scores[:2] # XXX

            # TODO base criteria on stdv of scores? 

            # Cut out obviously terrible options
            
            # TODO deal with case where makes infeasible
            # Could track which sites have been constrained in this manner and skip when they have.

            assert best_score>=0

            always_allow_second_best_score=True
            tolerable_score=np.inf 
            speedy=True
            if speedy:
                #tolerable_score=(best_score*1.5+50)  # Best when sigma is 1 from brief tests    
                # 
                # 
                #    
                use_second_best_score=False
                if use_second_best_score:
                    tolerable_score=second_best_score*2  # TODO try taking max of this value based off second best and one based off best
                else:
                    #tolerable_score=(best_score*20+1e3/3*n)
                    tolerable_score=(best_score*12+1e3/3*n)
                if always_allow_second_best_score:
                    tolerable_score=max(second_best_score,tolerable_score)

                #tolerable_score=np.inf
                #tolerable_score=best_score 
            debug_force_no_change=False
            if debug_force_no_change:
                tolerable_score=0
            '''
        else:
            assert False
            group_choices = [disordered_connection]
            perm_scores = [0]

        """

        nonlocal num_allowed_connections
        nonlocal num_forbidden_connections


        site_names = "|".join([f"{ch.resnum}.{ch.name}" for ch in disordered_connection[0].atom_chunks])


        if ordered_connection_option.hydrogen_tag!="":
            extra_tag = "Htag["+ordered_connection_option.hydrogen_tag+"]_"
        site_names+=extra_tag

        
        # allowed_connections = []
        # for p, (group,group_vars, score,code) in enumerate(zip(group_choices,group_choice_var_sets,perm_scores,group_choice_codes)):
            

        #     # So possible solution is guaranteed, always allow connections of original structure.  
        #     connection_altlocs = [
        #         {ch.altloc for ch in connection.atom_chunks} for connection in group
        #     ]
        #     is_no_swap_perm = all(len(altlocs_present) == 1 for altlocs_present in connection_altlocs)

        #     allowed = (score <= tolerable_score) or  is_no_swap_perm



        #     num_allowed_connections+=allowed 
        #     num_forbidden_connections+= not allowed



        #     tag=code
        #     if not allowed: 
        #         #pass
        #         lp_problem += (
        #             lpSum(group_vars)<=len(group_vars)-2,   # If all but one of the group vars is active, then the other must be active too. Hence we can do -2 not -1.
        #             f"FORBID{constraint_type.value}_{site_names}_{code}"
        #         )
        #     else:
        #         for ordered_connection_option,var_active in zip(group,group_vars):
        #             if ordered_connection_option not in allowed_connections:
        #                 allowed_connections.append(ordered_connection_option)
        #                 #allowed_connections.append((ordered_connection_option,var_active))

        #     # TODO if len altlocs is 2, then remove both connection options from allowed_connections.





        #     # This makes things slower...
        #     '''
        #     if allowed:
        #         assert len(permutation_vars)==len(altlocs)
        #         lp_problem += (
        #             lpSum(permutation_vars) <= len(permutation_vars)*permutation_vars[0],   # all or none are active.
        #             f"Permutation{constraint_type.value}{p}_{tag}"
        #             )
        #     '''
                    
            #TODO CRITICAL Variables for when hydrogen is involved can be simplified, especially for angles involving hydrogen.

        # TODO!!!! if connection options don't appear in any non-forbidden group, can just write:
        # for to_altloc in all_altlocs:
        #     assignment_vars = [site_var_dict[VariableID.Atom(ch)][ch.get_altloc()][to_altloc] for ch in ordered_connection_option.atom_chunks]
        #     lp_problem += (
        #         lpSum(assignment_vars) <=  len(assignment_vars)-1,   
        #         f"FORBID{constraint_type.value}_{tag}>>{to_altloc}"
        #         )
        # And also remove corresponding var_active ordered connection variable.
        # This will make things clearer, if not faster. 
        
        worst_no_change_score=0
        for altlocs_key in connection_var_dict:
            ordered_connection_option= connection_dict[altlocs_key]
            no_change_connection = len({ch.altloc for ch in ordered_connection_option.atom_chunks})==1
            if no_change_connection:
                worst_no_change_score=max(worst_no_change_score,ordered_connection_option.ts_distance)
        
        local_score_tolerate_threshold=10*worst_no_change_score
        always_tolerate_score_threshold = max(local_score_tolerate_threshold,global_score_tolerate_threshold) #worst_no_change_score*10+1e4

        for altlocs_key in connection_var_dict:
            #for ordered_connection_option,var_active in disordered_connection:
            ordered_connection_option,var_active = connection_dict[altlocs_key], connection_var_dict[altlocs_key]



            no_change_connection = len({ch.altloc for ch in ordered_connection_option.atom_chunks})==1
            allowed =  (not ordered_connection_option.forbidden) \
                or ordered_connection_option.ts_distance<=always_tolerate_score_threshold
            
            chance_allow_anyway=0
            if not allowed and chance_allow_anyway>0:
                if np.random.rand()<chance_allow_anyway:
                    allowed=True
                
            from_ordered_atoms = "|".join([f"{ch.resnum}.{ch.name}_{ch.altloc}" for ch in ordered_connection_option.atom_chunks])
            tag=from_ordered_atoms
            extra_tag=""
            if ordered_connection_option.hydrogen_tag!="":
                extra_tag = "Htag["+ordered_connection_option.hydrogen_tag+"]_"
            tag+=extra_tag

                
            # Connections contains ordered atoms from any and likely multiplke altlocs that 
            # *are to be assigned to the same altlocs*
            assignment_options:dict[str,list[LpVariable]]={}
            for to_altloc in all_altlocs:
                assignment_options[to_altloc]=[]
                for ch in ordered_connection_option.atom_chunks:
                    from_altloc = ch.get_altloc()
                    site = VariableID.Atom(ch)

                    if to_altloc not in site_var_dict[site][from_altloc]:
                        pass
                        #  print(f"Warning: to_altloc {to_altloc} not found in {site_var_dict[site][from_altloc]}, skipping {disordered_connection[0].get_disordered_connection_id()}")
                        # return

                    assignment_options[to_altloc].append(site_var_dict[site][from_altloc][to_altloc])
            # if variable is inactive, cannot have all atoms assigned to the same altloc.
            # Note that since every connection option is looped through, this also means 
            # that if variable is active, all atoms will be assigned to the same altloc.
            
            if allowed:
                distance_vars.append(ordered_connection_option.ts_distance*var_active)
            num_allowed_connections+=allowed 
            num_forbidden_connections+=not allowed 
            
            # Constraint will be handled by parent atom of hydrogen

            #disordered_connection_vars.append(var_active)

            for to_altloc, assignment_vars in assignment_options.items():
                #swaps = '|'.join([f"{''.join(assignment.name.split('_')[1:])}" for assignment in assignment_vars])
                num_assignments = len(ordered_connection_option.atom_chunks)
                assert len(assignment_vars) == num_assignments
                
                if allowed:
                    lp_problem += (
                        lpSum(assignment_vars) <=  num_assignments-1+var_active,   # Active if all assignment vars active.
                        f"ALLOW{constraint_type.value}_{tag}>>{to_altloc}"
                    )
                else: 
                    lp_problem += (
                        lpSum(assignment_vars) <=  num_assignments-1,   
                        f"FORBID{constraint_type.value}_{tag}>>{to_altloc}"
                    )
    
    worst_connection_before_swap=None
    worst_global_no_change_score=0
    for connection_id, ordered_connection_choices in disordered_connections.items():
        for c in ordered_connection_choices:
            preswap_connection=False
            if len({ch.altloc for ch in c.atom_chunks})==1:
                preswap_connection=True
            if preswap_connection and c.ts_distance> worst_global_no_change_score:
                worst_global_no_change_score = c.ts_distance
                worst_connection_before_swap = c

        
    # TODO make it the 25th percentile or something.
    global_score_tolerate_threshold=worst_global_no_change_score/100
    
    for i, (connection_id, ordered_connection_choices) in enumerate(disordered_connections.items()):
        if i % 250 == 0:
            print(f"Adding constraints {i}/{len(disordered_connections)}")
        constraint_type = VariableKind[connection_id.split('_')[0]]  #XXX ?????
        add_constraints_from_disordered_connection(constraint_type,ordered_connection_choices,global_score_tolerate_threshold=global_score_tolerate_threshold)
    print(f"Num allowed connections: {num_allowed_connections} | num forbidden connections: {num_forbidden_connections}")

    # if  force_sulfur_bridge_swap_solutions \
    #     and [ch.element for ch in connection.atom_chunks]==["S","S"]:
    #     forced_swap_solutions.append(
    #         var_dict[atom_a]["flipped"]==var_dict[atom_b]["flipped"]==1
    #     )

    if  force_sulfur_bridge_swap_solutions:
        raise Exception("force swap solutions unimplemented")







    badness = lpSum([dist for dist in distance_vars])
    lp_problem += (
    badness,
    "badness",)
    initial_badness = badness.value()
    print(f"Initial badness: {initial_badness}")
        
    
    # TODO: Suspect we want the two badness to be close to equal. Because if they aren't similar badness, we wouldn't expec to see both?

    
    ######################
    log_file = f"{out_dir}/xLO-HighLevelLog_{out_handle}.log"
    if not os.path.exists(log_file):
        with open(log_file,'w') as f:
            f.write(f"{out_handle} altloc optimizer log\n")
    def log(line:str):
        with open(log_file, 'a') as f:
            f.write(line+"\n")
            

    swaps_file =  swaps_file_path(out_dir,out_handle)
    def update_swaps_file(distances, site_assignment_arrays,record_notable_improvements_threshold=None): # record_notable_improvements_threshold: fractional improvement required to record separately
    # Create json file that lists all the site *changes* that are required to meet the solution. 
        out_dict = {"target": out_handle,"initial badness":initial_badness,"solutions":{}}
        if record_notable_improvements_threshold is not None:
            assert 0 <= record_notable_improvements_threshold < 1
            sep_dict= deepcopy(out_dict)
            best_improvement = 0
        verbose=False
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
                        # No change needed
                        continue
                    if verbose:
                        # TODO Flag only when swaps have **changed**, as did in 2 altloc version.
                        print(f"Flagging assignment of {site} {from_altloc} to {to_altloc}")
                    if site_key not in moves:
                        moves[site_key] = {}
                    moves[site_key][from_altloc] = to_altloc
        with open(swaps_file,'w') as f: 
            json.dump(out_dict,f,indent=4)

        if len(sep_dict["solutions"])>0:
            separate_record_file = f"{out_dir}/xLO-Diff_{f'{best_improvement*100:.2f}'}_{out_handle}.json"
            with open(separate_record_file,'w') as f2: 
                json.dump(sep_dict,f2,indent=4)


    assert num_solutions >= len(forced_swap_solutions) 

    site_assignment_arrays:list[dict[VariableID,dict[str,str]]]=[]
    distances=[]
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
            print(f"Solving {l}th best solution")
        print(f"-----------------------------------------------------")
        print()

        lp_problem.writeLP(f"{out_dir}/xLO-LP_{out_handle}.lp")


        class Solver(Enum):
            COIN="COIN"
            CPLX_PY="CPLX_PY"
            CPLX_CMD="CPLX_CMD"

        timeLimit=None
        if max_mins_start is not None:
            extra_time_per_loop=60*mins_extra_per_loop
            timeLimitStart=60*max_mins_start
            timeLimit=timeLimitStart+l*extra_time_per_loop
        elif mins_extra_per_loop != 0:
            print("Warning: extra time per loop is not 0, but there is no time limit")
        #timeLimit=None
        #threads=None
        #threads=26
        threads=10
        logPath=out_dir+"xLO-Log"+out_handle+".log"
        #logPath=None
        pulp_solver = Solver.CPLX_PY # https://stackoverflow.com/questions/10035541/what-causes-a-python-segmentation-fault
        #pulp_solver = Solver.COIN
        warmStart=True
        #gapRel=0.0003
        #gapRel=0.001
        
        #gapRel=None

        with open(f"{out_dir}/xLO-Initial{out_handle}.txt",'w') as f:
            f.write(f"Status: {LpStatus[lp_problem.status]}\n")

            for v in lp_problem.variables():
                try:
                    if v.value() > 0:
                        f.write(f"{v.name} = {v.value()}\n")
                except:
                    raise Exception(v,v.value())
            f.write(f"Total distance = {value(lp_problem.objective)}")


        solver_class=PULP_CBC_CMD
        solver_options=[]
        if pulp_solver == Solver.COIN: 
            solver_class = PULP_CBC_CMD
        elif pulp_solver == Solver.CPLX_CMD:   # cant get working
            assert False
            pulp_solver = CPLEX_CMD
        # CPLEX parameters: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/tutorials/InteractiveOptimizer/settingParams.html
        # CPLEX status: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.10.0/ilog.odms.cplex.help/refcallablelibrary/macros/Solution_status_codes.html
            solver_options.append("set parallel -1")
            #path='~/ibm/ILOG/CPLEX_STUDIO2211/cplex'
        elif pulp_solver == Solver.CPLX_PY:
            solver_class = CPLEX_PY
        else:
            raise Exception("not implemented")
        #solver = solver_class(timeLimit=timeLimit,threads=threads,logPath=logPath,warmStart=warmStart,path=path)
        solver = solver_class(timeLimit=timeLimit,threads=threads,warmStart=warmStart,logPath=logPath,gapRel=gapRel)
        lp_problem.solve(solver)
        print(solver)
        # if pulp_solver == Solver.COIN:
        #     # solver = PULP_CBC_CMD(threads=NTHREADS)
        #     # lp_problem.solve(solver=solver)
        #     maxNodes=None
        #     #pulpTestAll()
        #     #asdds
        #     solver = PULP_CBC_CMD(logPath=,threads=threads,timeLimit=timeLimit,maxNodes=maxNodes)
        # elif pulp_solver == Solver.COIN:
        #     solver = CPLEX_PY(timeLimit=timeLimit,threads=threads)
        # lp_problem.solve(solver)

        
        def get_status(verbose=False):
            print("Status:", LpStatus[lp_problem.status])
            

            if verbose:
                for v in lp_problem.variables():
                    if v.value() > 0:
                        print(v.name, "=", v.value())

            print(f"Target: {out_handle}")
            total_distance = value(lp_problem.objective)
            diff=total_distance/initial_badness-1
            print(f"Total distance = {total_distance} ({100*(diff):.3f}%)")
            log(f"{100*(diff):.3f}%")
            #plt.scatter()
        get_status(verbose=False)




        # dry = None
        # for val in connections.values():
        #     for v in val.values():
        #         dry=v.calculated_dry()
        #         break
        # tag = "_dry" if dry else ""
        with open(f"{out_dir}/xLO-Out{out_handle}.txt",'w') as f:
            f.write(f"Status: {LpStatus[lp_problem.status]}\n")

            for v in lp_problem.variables():
                try:
                    if v.value() > 0:
                        f.write(f"{v.name} = {v.value()}\n")
                except:
                    raise Exception(v,v.value())
            f.write(f"Total distance = {value(lp_problem.objective)}")


        if lp_problem.sol_status==LpStatusInfeasible:
            assert l>0, "Solution was infeasible!"

            print()
            print(f"WARNING: Finding solution {l} was infeasible! Ending solution search")
            break

        with open(f"{out_dir}/xLO-ActiveConnections{out_handle}.txt",'w') as f:
            f.write("name, sigma, cost\n")
            vals = [(constraint,var) for constraint,var in constraint_var_dict.values() if var.value()>0]
            vals.sort(key=lambda x: x[0].z_score,reverse=True)
            #vals.sort(key=lambda x: x[0].ts_distance)
            for constraint, var in vals:
                f.write(f"{var.name} {constraint.z_score:.2e} {constraint.ts_distance:.2e}\n")
        
        site_assignments:dict[VariableID,dict[str,dict[str,int]]] = {}
        site_assignment_arrays.append(site_assignments)



        # Determine which atom has been assigned where.
        for site in site_var_dict:
            site_assignments[site]={}
            for from_altloc in site_var_dict[site]:
                if lp_problem.sol_status==LpStatusInfeasible:
                    site_assignments[site][from_altloc]=from_altloc 
                    break 
                to_altloc_found = False
                for to_altloc in site_var_dict[site][from_altloc]:
                    if round(site_var_dict[site][from_altloc][to_altloc].value())==1:  # For some reason CPLEX outptuts values like 1.0000000000094025 sometimes.
                        assert not to_altloc_found
                        to_altloc_found=True
                        site_assignments[site][from_altloc]=to_altloc
                assert to_altloc_found, (site, from_altloc)

            # if site.name==str(DisorderedTag(10,"CA")) or site.name==str(DisorderedTag(2,"CA")) :
            #     for from_altloc in site_var_dict[site]:
            #         for to_altloc in site_var_dict[site][from_altloc]:

            #             print(from_altloc,to_altloc)
            #             print(round(site_var_dict[site][from_altloc][to_altloc].value()))
            #             print(site_var_dict[site][from_altloc][to_altloc].value())
            #             print(site_var_dict[site][from_altloc][to_altloc])


        # print("KEYS")
        # for key in nonzero_variables:
        #     print(key)
        distances.append(value(lp_problem.objective))
        
        update_swaps_file(distances,site_assignment_arrays,record_notable_improvements_threshold=0.03)
      
       
        # solution = [var.value() for var in lp_problem.variables()] 
        # lp_variables = lp_problem.variables()
        # lp_problem += (  
        #         pulp.lpSum(solution[i]*lp_variables[i] for i in range(len(solution))) <= len(solution)-1,
        #         f"forceNextBest_{l}"
        #     )


        flipped_flip_variables = []
        flip_variables=[]
        #TODO Critical - for some reason the 'next-best' solutions can be better. Maybe the gap tolerance is not 0?
        for chunk in chunk_sites.values():
            site = VariableID.Atom(chunk)
            if not site_being_considered(site):
                continue
            if protein_sites and site.is_water: # not interested in different solutions for water. 
                continue # This means water atoms may or may not swap for the single solution where no protein atoms swap.
            from_altloc = chunk.get_altloc()
            for to_altloc, var in site_var_dict[site][from_altloc].items():
                val = var.value()
                flip_variables.append(var)
                if val == 1:
                    #if (1-var) not in flipped_flip_variables:
                    flipped_flip_variables.append(1-var)
        if len(flipped_flip_variables) == 0:
            # Require one variable to be flipped
            lp_problem += pulp.lpSum(flip_variables) >= 1, f"force_swaps_loop_{l}"
        else:
            # require at least one flip to be different
            lp_problem += pulp.lpSum(flipped_flip_variables) >= 1, f"force_next_best_solution_{l}"
        
    del lp_problem
    del solver
    gc.collect()
        ##################

    assert swaps_file==swaps_file_path(out_dir,out_handle) # XXX
    return swaps_file

def swaps_file_path(out_dir,out_handle):
    return f"{out_dir}/xLO-toFlip_{out_handle}.json"

if __name__=="__main__":
    handle = sys.argv[1]
    out_dir = f"{os.path.abspath(os.getcwd())}/output/{handle}/"
    os.makedirs(out_dir,exist_ok=True)
    pdbPath = f"{UntangleFunctions.pdb_data_dir()}/{handle}.pdb"
    finest_chunks,disordered_connections= MTSP_Solver(pdbPath).calculate_paths(
        atoms_only=True,
    )
    solve(finest_chunks,disordered_connections,handle)


# %%
