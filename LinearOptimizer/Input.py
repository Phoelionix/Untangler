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

# Purpose: Get geometric badness for all possible "groups" of atoms that 
# we have geometric measures for, to pass to LinearOptimizer.Solver.


# TODO  
# - Change to proper chi-squared (missing denominator presently). 
# - implement all the other wonderful measures we have available. -- Both Input and Solver
# - Clean up chunk ID system. -- Both Input and Solver
# - Consider turning sigmas back off? Don't actually make much sense if following unrestrained refinement. Should divide by ideal instead.
# - If only some atoms are allowed, need to include all restraints that involve the atoms, even those that do not involve allowed atoms
# - # Swaps can create nonbond issues that are not recorded due to not being present in geo file?
 

from Bio.PDB import PDBParser,Structure
from Bio.PDB.Atom import Atom,DisorderedAtom
from Bio.PDB.Residue import Residue # Note we don't want the disorderedresidue here, it refers to different residue types for same res seq num.
from Bio.PDB.StructureBuilder import StructureBuilder
from pulp import *
import UntangleFunctions 
from UntangleFunctions import two_key_read
from LinearOptimizer.Swapper import Swapper
import os
import random
import numpy as np
import itertools
from multiprocessing import Pool
import scipy.stats as st
from statistics import NormalDist
from LinearOptimizer.VariableID import *
from LinearOptimizer.Tag import *
from LinearOptimizer.mon_lib_read import read_vdw_radii,read_lj_parameters,get_mon_lib_names
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler
from Untwist import untwist

#ALTERNATE_POSITIONS_FACTOR=1
ALTERNATE_POSITIONS_FACTOR=1
ALTERNATE_POSITIONS_Z_SCORE_FACTOR=1
#HYDROGEN_RESTRAINTS=False # If False, hydrogen restraints will be ignored.
NEVER_FORBID_HYDROGEN_GEOMETRIES=True
hydrogen_restraint_scale=1 # scales cost of hydrogen restraints. E.g. angles for serine OG - CB - H will prevent a swap that improves OG - CB bond length, even though after refine the hydrogens will rearrange with no issues.  
APPLY_TENSION_MOD=False
TURN_OFF_MIN_SIGMAS=False
ALLOW_OUTLIERS_FROM_POTENTIAL_OVERRIDES=False # Will not do min sigmas with  outliers identified as potentially genuine in {model_handle}_potential_overrides.txt.
#weight_mod_for_allowed_outliers=1e-1
weight_mod_for_allowed_outliers=1
#high_tension_penalty=1.75  # Penalty for current connections (i.e. geometries? groups?) that correspond to a high tension (above min_tension_where_anything_goes). TODO should just be a penalty for not having any geometries involving the site change.
#high_tension_penalty=10  # Penalty for current connections (i.e. geometries? groups?) that correspond to a high tension (above min_tension_where_anything_goes). TODO should just be a penalty for not having any geometries involving the site change.
high_tension_penalty=1  # Penalty factor for current connections (i.e. geometries? groups?) that correspond to a high tension (above min_tension_where_anything_goes). TODO should just be a penalty for not having any geometries involving the site change.

#TODO could try reverting model back to preswap model.

assert high_tension_penalty>=1



def is_atom(atom:Atom,resnum,name,altloc):
        if (atom.name==name and atom.get_altloc()==altloc and OrderedAtomLookup.atom_res_seq_num(atom)==resnum):
            return True
        return False
def atoms_in_LO_variable_string(variable:str,atoms:list[Atom]): # e.g. "30.C_B|31.CA_B" or "Nonbond_30.CA_A|31.N_A"
    if not variable.split("_")[0].isnumeric():
        variable = variable[variable.find("_")+1:]
    atom_strings = variable.split("|")
    
    for a in atoms:
        for atom_string in atom_strings:
            resnum, name_altloc = atom_string.split('.')
            resnum = int(resnum)
            name, altloc = name_altloc.split("_")  
            if is_atom(a,resnum,name,altloc):
                break
        else: # did not find a match
            return False
    return True

# Want to calculate wE for chunks of atoms and residues and side chains and main chains 
# we then get the optimizer to minimize when choosing 
# E.g. calculate chunk over two halves of residue
# So if route includes those halves, the cost includes those.
class Chunk():
    def __init__(self,depth_tag,site_num,altloc,atoms:list[Atom],start_resnum,end_resnum,start_resname,end_resname):
        self.altloc=altloc
        self.site_num_at_depth=site_num # Chunks at same disordered site have same site num
        self.atoms_dict:dict[int,list[Atom]]={}
        for atom in atoms:
            res_num = OrderedAtomLookup.atom_res_seq_num(atom)
            if res_num not in atoms:
                self.atoms_dict[res_num] = []
            self.atoms_dict[res_num].append(atom)
        self.start_resname=start_resname
        self.end_resname=end_resname
        self.start_resnum = start_resnum
        self.end_resnum = end_resnum
        self.depth_tag=depth_tag
    def get_site_num(self):
        return self.site_num_at_depth
    def get_altloc(self):
        return self.altloc
    def get_disordered_tag(self):
        pass
        # Making it clear this code isn't run because I made this so confusing. But might be useful idea e.g. for van der waals
        #return  f"{self.start_resnum}.{self.start_resname}x{self.end_resnum}.{self.end_resname}"
    def unique_id(self):
        pass
        #return f"{self.depth_tag}&{self.get_site_num()}.{self.get_altloc()}&{self.get_disordered_tag()}"
    
    def get_atoms_dict(self):
      return self.atoms_dict


class OrderedResidue(Chunk):
    def __init__(self,altloc,resnum,referenceResidue,atoms,site_num_override=None):
        self.resnum=resnum
        self._res=referenceResidue
        site_num = self.get_resnum() # Index of disordered 
        if site_num_override is not None:
            site_num =site_num_override
        super().__init__("RESI",site_num,altloc,atoms,
                         self.get_resnum(),self.get_resnum(),
                         self.get_resname(),self.get_resname())
    def get_resnum(self):
        return self.resnum
    def get_resname(self):
        return self._res.get_resname()
    def get_id(self):
        return self._res.get_id()
    def get_atoms(self):
        return self.atoms_dict[self.get_resname()]
    def get_disordered_tag(self):
        assert False 
        return  f"{self.get_resnum()}.{self.get_resname()}"
    

# Ugh, never do inheritance. TODO refactor to composition.
class AtomChunk(OrderedResidue): 
    # Just an atom c:
    def __init__(self,site_num,altloc,resnum,referenceResidue,atom:Atom,is_water:bool,constraints_handler,alt_pos_options:List[tuple[int,'Any']]):
        self.name=atom.get_name()
        self.is_water = is_water
        self.element = atom.element
        self.coord = atom.get_coord()
        self.alt_pos_options:dict[int,Atom] = alt_pos_options
        super().__init__(altloc,resnum,referenceResidue,[atom],site_num)
        self.generate_uninitialized_constraints(constraints_handler)
        self.depth_tag="ATOM"
    def generate_uninitialized_constraints(self,constraints_handler:ConstraintsHandler):
        self.constraints_holder:ConstraintsHandler = \
            constraints_handler.get_constraints(self.get_disordered_tag(),no_constraints_ok=self.is_water)
    def has_alternate_coords(self):
        return len(self.alt_pos_options)>0
    def get_coord(self):
        return self.coord
    def get_disordered_tag(self):
        return DisorderedTag(self.get_resnum(),self.name)
    def get_ordered_tag(self):
        return OrderedTag(self.get_resnum(),self.name,self.altloc)
    def unique_id(self):
        return f"{self.depth_tag}&{self.get_site_num()}.{self.get_altloc()}&{self.get_disordered_tag()}"


class LP_Input:
    #MODE="LOW_TOL" # "NONBOND_RESTRICTIONS" #"LOW_TOL" #"HIGH_TOL" #"NO_RESTRICTIONS" # PHENIX REFMAC
    #MODE= "NO_RESTRICTIONS"
    MODE= "HIGH_TOL"

    if MODE=="NO_RESTRICTIONS":
        max_sigmas=min_sigmas_where_anything_goes=min_tension_where_anything_goes={}

    elif MODE=="NONBOND_RESTRICTIONS":
        max_sigmas={
            ConstraintsHandler.NonbondConstraint:3,
        }
        min_sigmas_where_anything_goes={
            ConstraintsHandler.NonbondConstraint:2,
        }
        min_tension_where_anything_goes={}
    elif MODE=="HIGH_TOL":
        max_sigmas={
            ConstraintsHandler.BondConstraint:10,
            ConstraintsHandler.AngleConstraint:4,
        }    
        min_sigmas_where_anything_goes={
            ConstraintsHandler.BondConstraint:2,
            ConstraintsHandler.NonbondConstraint:2,
        } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:5,
            ConstraintsHandler.AngleConstraint:5,
            ConstraintsHandler.NonbondConstraint:5,
            ConstraintsHandler.ClashConstraint:5,
            # ConstraintsHandler.BondConstraint:7,
            # ConstraintsHandler.AngleConstraint:7,
            # ConstraintsHandler.NonbondConstraint:7,
            # ConstraintsHandler.ClashConstraint:7,
        } 
    elif MODE=="MED_TOL":
        max_sigmas={
            ConstraintsHandler.BondConstraint:8,
            ConstraintsHandler.AngleConstraint:3,
        }    
        min_sigmas_where_anything_goes={
            #ConstraintsHandler.BondConstraint:99,
            #ConstraintsHandler.AngleConstraint:99,
            #ConstraintsHandler.NonbondConstraint:2,
        } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:8,
            ConstraintsHandler.AngleConstraint:4,
            #ConstraintsHandler.BondConstraint:8,
            #ConstraintsHandler.AngleConstraint:8,
            ConstraintsHandler.NonbondConstraint:8,
            ConstraintsHandler.ClashConstraint:8,
        } 
    elif MODE=="LOW_TOL":
        max_sigmas={
            ConstraintsHandler.BondConstraint:8,
            ConstraintsHandler.AngleConstraint:2.5,
        }    
        min_sigmas_where_anything_goes={
            #ConstraintsHandler.BondConstraint:99,
            #ConstraintsHandler.AngleConstraint:99,
            #ConstraintsHandler.NonbondConstraint:2,
        } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:8,
            ConstraintsHandler.AngleConstraint:4,
            #ConstraintsHandler.BondConstraint:8,
            #ConstraintsHandler.AngleConstraint:8,
            ConstraintsHandler.NonbondConstraint:8,
            ConstraintsHandler.ClashConstraint:8,
        } 
    elif MODE=="TENSIONS_TOL":
        max_sigmas={
            ConstraintsHandler.BondConstraint:2.5,
            ConstraintsHandler.AngleConstraint:2.5,
            ConstraintsHandler.NonbondConstraint:2.5,
        }    
        min_sigmas_where_anything_goes={
        } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:2.5,
            ConstraintsHandler.AngleConstraint:2.5,
            ConstraintsHandler.NonbondConstraint:2.5,
            ConstraintsHandler.ClashConstraint:2.5,
        } 
    elif MODE=="PHENIX":
        max_sigmas={
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            #ConstraintsHandler.BondConstraint:6,
            #ConstraintsHandler.AngleConstraint:3,
            ConstraintsHandler.BondConstraint:4,
            ConstraintsHandler.AngleConstraint:3,
            # ConstraintsHandler.BondConstraint:2,
            # ConstraintsHandler.AngleConstraint:2,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        min_sigmas_where_anything_goes={
            ConstraintsHandler.BondConstraint:2.5,
            ConstraintsHandler.AngleConstraint:2.5,
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        # min_tension_where_anything_goes={
        #     ConstraintsHandler.BondConstraint:4,
        #     ConstraintsHandler.AngleConstraint:5,
        #     ConstraintsHandler.NonbondConstraint:4,
        #     ConstraintsHandler.ClashConstraint:4,
        # } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:6,
            ConstraintsHandler.AngleConstraint:6,
            ConstraintsHandler.NonbondConstraint:6,
            ConstraintsHandler.ClashConstraint:6,
            # ConstraintsHandler.BondConstraint:7,
            # ConstraintsHandler.AngleConstraint:7,
            # ConstraintsHandler.NonbondConstraint:7,
            # ConstraintsHandler.ClashConstraint:7,
        } 
    elif MODE == "REFMAC":
        max_sigmas={
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            # ConstraintsHandler.BondConstraint:2.5,
            # ConstraintsHandler.AngleConstraint:2.5,
            ConstraintsHandler.BondConstraint:4,
            ConstraintsHandler.AngleConstraint:4,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        min_sigmas_where_anything_goes={
            # ConstraintsHandler.BondConstraint:2.5,
            # ConstraintsHandler.AngleConstraint:2.5,
            ConstraintsHandler.BondConstraint:2,
             ConstraintsHandler.AngleConstraint:3,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        # min_tension_where_anything_goes={
        #     ConstraintsHandler.BondConstraint:4,
        #     ConstraintsHandler.AngleConstraint:5,
        #     ConstraintsHandler.NonbondConstraint:4,
        #     ConstraintsHandler.ClashConstraint:4,
        # } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:6,
            ConstraintsHandler.AngleConstraint:6,
            ConstraintsHandler.NonbondConstraint:6,
            ConstraintsHandler.ClashConstraint:6,
            # ConstraintsHandler.BondConstraint:5,
            # ConstraintsHandler.AngleConstraint:5,
            # ConstraintsHandler.NonbondConstraint:5,
            # ConstraintsHandler.ClashConstraint:5,
        } 
    else:
        raise Exception(f"Invalid MODE {MODE}")
    ### This commented out code computes wE for combinations of larger groups of atoms.  
    # Could be useful in future as a quick coarse step... but would probably be better 
    # to build up from the atom-site approach, coding the wE measure directly.   
    ### 

    # class ChunkConnection():
    #     [ deleted ]
    #         _,self.ts_distance,_ = UntangleFunctions.assess_geometry_wE(connection_structure_save_path,tmp_out_folder_path,phenixgeometry_only=quick_wE) 


    @staticmethod
    def make_hydrogen_tag(hydrogen_names):
        return "_"+''.join(hydrogen_names) if hydrogen_names is not None else ""
    
    # Connection between ordered atoms
    class Geomection(): # TODO really need to have a disordered connection class to have some of these properties (e.g. max site tension, outlier_ok)

        def __init__(self, atom_chunks:list[AtomChunk],ts_distance,position_option_indices:list[int],connection_type,hydrogen_names,ideal,z_score,max_site_tension,outlier_ok):
            assert hydrogen_names==None or len(hydrogen_names)==len(atom_chunks)
            self.atom_chunks = atom_chunks
            self.from_altlocs=[a.get_altloc() for a in atom_chunks]
            self.res_nums=[a.get_resnum() for a in atom_chunks]
            self.atom_names=[a.name for a in atom_chunks]
            self.connection_type=connection_type
            self.ts_distance=ts_distance  # NOTE As in the travelling salesman problem sense
            self.position_option_indices=position_option_indices
            self.hydrogen_tag = LP_Input.make_hydrogen_tag(hydrogen_names)
            self.poschange_tag=""
            self.hydrogen_name_set=set([])
            self.ideal=ideal
            self.z_score=z_score # i.e. sigma
            ###
            self.max_site_tension=max_site_tension
            self.outlier_ok=outlier_ok
            ###
            #TODO move this to LinearOptimizer.Solver.py!!! 
            self.forbidden= (connection_type in LP_Input.max_sigmas) and (self.z_score > LP_Input.max_sigmas[self.connection_type]) 
                
            if hydrogen_names is not None:
                self.hydrogen_name_set = set(hydrogen_names)
                if NEVER_FORBID_HYDROGEN_GEOMETRIES:
                    self.forbidden=False
            if self.involves_position_changes():
                self.poschange_tag = '_'+'.'.join([str(i) for i in self.position_option_indices])
        def involves_position_changes(self):
            return any([i != 0 for i in self.position_option_indices])
        def single_altloc(self):
            return len({ch.altloc for ch in self.atom_chunks})==1
        # Whether geomection is active in preswap structure
        def original(self):
            return self.single_altloc() and not self.involves_position_changes()
        def get_disordered_connection_id(self):
            kind = ConstraintsHandler.Constraint.kind(self.connection_type)
            return f"{kind}{self.hydrogen_tag}_{'_'.join([str(a_chunk.get_disordered_tag()) for a_chunk in self.atom_chunks])}"

        # TODO TEMPORARY XXX
        @staticmethod
        def construct_disordered_connection_id(connection_type,disordered_tags:list[DisorderedTag],hydrogen_names=None):
            kind = ConstraintsHandler.Constraint.kind(connection_type)
            hydrogen_tag = LP_Input.make_hydrogen_tag(hydrogen_names)
            return f"{kind}{hydrogen_tag}_{'_'.join([str(dtag) for dtag in disordered_tags])}"


    def __init__(self,pdb_file_path:str,restrained_refine_pdb_file_path:str,APPLY_TENSION_MOD:dict[DisorderedTag,float], symmetries, align_uncertainty=False,ignore_waters=False,altloc_subset=None,resnums=None,resnames=None,
                 alternate_atoms=[]): 
        # TODO when subset size > 2, employ fragmentation/partitioning.
         
        # Note if we ignore waters then we aren't considering nonbond clashes between macromolecule and water.
        original_structure = PDBParser().get_structure("struct",pdb_file_path)
        # if align_uncertainty:
        #     self.align_uncertainty(original_structure)
        
        # debug_quick=False
        # if debug_quick:
        #     resnums = range(64)


        # Disabling because phenix seems to be deleting remark 290 :C
        '''
        self.symmetries = UntangleFunctions.parse_symmetries_from_pdb(pdb_file_path)
        '''
        self.symmetries:list = symmetries

        self.ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=not ignore_waters,
                                                     altloc_subset=altloc_subset,
                                                     allowed_resnums=resnums,allowed_resnames=resnames,
                                                     alternate_atoms=alternate_atoms)   
        self.model_path=LP_Input.subset_model_path(pdb_file_path,altloc_subset)
        self.restrained_model_path=restrained_refine_pdb_file_path # TODO make optional arg default None. And if None, then in calculate_paths, set original clashes to [].
        self.APPLY_TENSION_MOD=APPLY_TENSION_MOD
        
    def align_uncertainty(self,structure:Structure.Structure):
        # in x-ray data and geom.
        # Experimental and doesn't help.
        for atom in structure.get_atoms():
            a,b = atom.__iter__()
            coord_diff = np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
            mean_coord = (a.get_coord()+b.get_coord())/2
            if coord_diff < 0.08: #TODO resolution based #TODO properly do the rotation and find minimum method
                a.coord = b.coord =  mean_coord
                atom.coord = mean_coord

    @staticmethod 
    def subset_model_path(pdb_file_path,altloc_subset:list[str]):
        tag=""
        if altloc_subset is not None:
            if len(altloc_subset)==1:
                tag = f"_conformer-{''.join(altloc_subset)}"
            else:
                tag = f"_subset-{''.join(altloc_subset)}"
        assert pdb_file_path[-4:]==".pdb",pdb_file_path
        return os.path.join(UntangleFunctions.separated_conformer_pdb_dir(), os.path.basename(pdb_file_path)[:-4]+tag+".pdb")
    @staticmethod
    def prepare_geom_files(base_model_path,all_altloc_subsets:list[str],num_threads=10,water_swaps=False,allowed_resnums=None,allowed_resnames=None,waters=True):
        if all_altloc_subsets is None:
            all_altloc_subsets=[None]
        if all_altloc_subsets==[None]:
            print("Considering full set of altlocs")
        elif (type(all_altloc_subsets)==str) or len(all_altloc_subsets[0])==1:
            raise Exception(f"all_altloc_subsets is {all_altloc_subsets}. But it should be either 1. a list of altloc subsets (strings with more than 1 character/altloc) to prepare files for, or 2. None, to consider all altlocs")

        original_structure = PDBParser().get_structure("struct",base_model_path)
        subset_model_paths=[]

        for altloc_subset in all_altloc_subsets:
            ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=waters,
                                                     altloc_subset=altloc_subset,
                                                     allowed_resnums=allowed_resnums,allowed_resnames=allowed_resnames)
                    
            subset_model = LP_Input.subset_model_path(base_model_path,altloc_subset)
            assert subset_model not in subset_model_paths
            subset_model_paths.append(subset_model)
            ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=base_model_path,out_path=subset_model)


        global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
        def pooled_method(i):
            LP_Input.prepare_geom_files_for_one_subset(subset_model_paths[i],water_swaps=water_swaps)

        with Pool(num_threads) as p:
            p.map(pooled_method,range(len(subset_model_paths)))
            

    ''' Assess conformations separately. TODO implement
    def geo_model_paths(pdb_file_path,altloc_subset:list[str]):
        assert pdb_file_path[-4:]==".pdb",pdb_file_path
        assert altloc_subset is not None
        files=[]
        for altloc in altloc_subset:
            tag = f"_conformer-{altloc}"
            files.append(os.path.join(
                UntangleFunctions.separated_conformer_pdb_dir(), os.path.basename(pdb_file_path)[:-4]+tag+".pdb"
            ))
        return files
            

    @staticmethod
    def prepare_geom_files(base_model_path,all_altloc_subsets,num_threads=10,water_swaps=False):

        original_structure = PDBParser().get_structure("struct",base_model_path)

        # Create subset pdb files
        for altloc_subset in all_altloc_subsets:
            ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=True,
                                                     altloc_subset=altloc_subset)
                    
            subset_model = LP_Input.subset_model_path(base_model_path,altloc_subset)
            ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=base_model_path,out_path=subset_model)

        altlocs = []
        for altloc_subset in all_altloc_subsets:
            for altloc in altloc_subset:
                if altloc not in altlocs:
                    altlocs.append(altloc)

        # Create geo files for each conformation
        geo_model_paths=[]
        for altloc in altlocs:
            ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=True,
                                                     altloc_subset=[altloc])
                    
            conformation = LP_Input.geo_model_paths(base_model_path,[altloc])
            assert len(conformation)==1 and type(conformation)==list
            conformation=conformation[0]
            assert conformation not in geo_model_paths
            geo_model_paths.append(conformation)
            ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=base_model_path,out_path=conformation)
        

        global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
        def pooled_method(i):
            LP_Input.prepare_geom_files_for_one_subset(subset_model_paths[i],water_swaps=water_swaps)

        with Pool(num_threads) as p:
            p.map(pooled_method,range(len(subset_model_paths)))
    '''


    @staticmethod
    def geo_log_out_folder():
        return UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+"StructureGeneration/HoltonOutputs/"
    
    @staticmethod
    def prepare_geom_files_for_one_subset(model_path,water_swaps=False,turn_off_cdl=False): # TODO CDL incorporation.
        UntangleFunctions.assess_geometry_wE(model_path,turn_off_cdl=turn_off_cdl) 

    @staticmethod 
    def water_swapped_handle(model_path):
        return UntangleFunctions.model_handle(model_path)+"_WaSw"
        
    def calculate_paths(self,scoring_function,quick_wE=False, dry_run=False,atoms_only=True,
                        clash_punish_thing=False,nonbonds=True,water_water_nonbond=None,
                        constraint_weights:dict[Type,float]=None, # dict of serial numbers with altloc labels. All with same label will be put in same conformation. The actual label of each conformation in the output file may be different. This only specifies which are to be in the same conformation.
                        force_solution_reference:dict[OrderedTag,str]=None,
                        suppress_worse=False, # Suppresses cost of ordered connections for a geometry when they are worse than the current worst ordered connection. This is to encourage improvements in other measures. 
                        improvement_factor=1, # improvements are weighted by this much. Might need to be implemented in LinearOptimizer.Solver
                        weight_for_range=False, # increases weight of restraints that currently have a wide range in z scores. 
                        )->tuple[list[Chunk],dict[str,list[Geomection]]]: #disorderedResidues:list[Residue]
        print("Calculating geometric costs for all possible connections between chunks of atoms (pairs for bonds, triplets for angles, etc.)")
        
        if constraint_weights is None:
            constraint_weights={}
        
        for key in [
            ConstraintsHandler.BondConstraint,
            ConstraintsHandler.AngleConstraint,
            ConstraintsHandler.NonbondConstraint,
            ConstraintsHandler.ClashConstraint,
            ConstraintsHandler.TwoAtomPenalty,
        ]:
            if key not in constraint_weights:
                constraint_weights[key]=1

        tmp_out_folder_path=UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+"LinearOptimizer/tmp_out/"
        os.makedirs(tmp_out_folder_path,exist_ok=True)
        
        assert atoms_only, "solving with wEs of chunks of residues not supported right now"

        chunk_sets = []
        connection_types=[]
        if not atoms_only:
            assert False
            orderedResidues: list[OrderedResidue]= [] # residues for each conformation
            self.dry_run=dry_run
            for res_num,referenceResidue in zip(self.ordered_atom_lookup.get_residue_nums(),self.ordered_atom_lookup.get_residue_sources()):
                assert not UntangleFunctions.res_is_water(referenceResidue)
                for altloc in self.ordered_atom_lookup.get_altlocs():
                    atoms = self.ordered_atom_lookup.select_atoms_by(res_nums=[res_num],altlocs=[altloc])
                    orderedResidues.append(OrderedResidue(altloc,res_num,referenceResidue,atoms))
            chunk_sets.append(orderedResidues)
            connection_types.append(LP_Input.ChunkConnection)


            #Triplets of residues
            tripletResidues: list[Chunk]= []
            n_size=3
            for i,start_resnum in enumerate(self.ordered_atom_lookup.get_residue_nums()[::n_size]):
                end_resnum=min(start_resnum+n_size-1,max(self.ordered_atom_lookup.get_residue_nums()))
                start_resname = self.ordered_atom_lookup.res_names[start_resnum]
                end_resname = self.ordered_atom_lookup.res_names[end_resnum]
                for altloc in self.ordered_atom_lookup.get_altlocs():
                    atoms = self.ordered_atom_lookup.select_atoms_by(res_nums=range(start_resnum,end_resnum+1),altlocs=[altloc])
                    triplet = Chunk("TRPL",int((i+2)/2),altloc,atoms,start_resnum,end_resnum,start_resname,end_resname)
                    tripletResidues.append(triplet)
            # half-halfs of residues
            chunk_sets.append(tripletResidues)
            connection_types.append(LP_Input.ChunkConnection)

            tripletAtoms = "TODO"
        
        nonbond_scores_path=nonbond_water_flipped_scores_path=None
        model_handle = UntangleFunctions.model_handle(self.model_path)


        if water_water_nonbond is None:
            water_water_nonbond = self.ordered_atom_lookup.waters_allowed
        if not self.ordered_atom_lookup.waters_allowed and (water_water_nonbond):
            assert False

        # Generate geo file
        nonbond_scores_files = []
        clashes=[]
        if nonbonds:

            # CLASH 26.6682 0.79 | A  62 AARG  HD2| S 128 AHOH  O 
            # Terrible code. XXX
            def get_water_clashes(handle,from_altloc_dict:bool)->list: # from_altloc_dict, with keys being the to_altlocs.
                clashes = []
                #nonbond_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_scorednonbond.txt"
                clashes_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_clashes.txt"
                with open(clashes_path) as clashes_file:
                    for line in clashes_file:
                        #pdb1 = f"{name}     ARES     A      {res_num}"
                        res_name = [entry.split()[-2][1:] for entry in line.split("|")[1:]] 
                        if not res_name.count("HOH")==1:
                            continue
                        res_num = [int(entry.split()[1]) for entry in line.split("|")[1:]] 
                        #  NOTE At the moment altlocs will all be the same. But hope to 
                        # one day read from file that does clashes between all altlocs!!!!
                        to_altloc = [entry.split()[-2][0] for entry in line.split("|")[1:]] 
                        name = [entry.strip().split()[-1] for entry in line.split("|")[1:]] 
                        from_altloc = [a for a in to_altloc]
                        water_idx = res_name.index("HOH") # Because we continue above if don't have exactly one HOH involved in clash.
                        from_altloc[water_idx] = from_altloc_dict[to_altloc[water_idx]]
                        badness = float(line.split()[1])
                        #new_line = f"CLASH   {badness} XXXXX XXXXX XXXXX X |  {name[0]}_{res_num[0]} {name[1]}_{res_num[1]}"
                        clashes.append((name, res_num, badness, from_altloc))
                return clashes

            def get_clashes(handle)->list: # from_altloc_dict, with keys being the to_altlocs.
                clashes = []
                #nonbond_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_scorednonbond.txt"
                clashes_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_clashes.txt"
                print(f"Reading clashes from {clashes_path}")
                assert os.path.exists(clashes_path), f"{clashes_path} does not exist." 
                assert os.path.isfile(clashes_path), f"{clashes_path} is a directory!" 
                with open(clashes_path) as clashes_file:
                    for line in clashes_file:
                        #pdb1 = f"{name}     ARES     A      {res_num}"
                        res_num = [int(entry.split()[1]) for entry in line.split("|")[1:]] 
                        #  NOTE Altlocs will all be the same. 
                        altlocs = [entry.split()[-2][0] for entry in line.split("|")[1:]] 
                        
                        name = [entry.strip().split()[-1] for entry in line.split("|")[1:]] 
                        badness = float(line.split()[1])
                        #new_line = f"CLASH   {badness} XXXXX XXXXX XXXXX X |  {name[0]}_{res_num[0]} {name[1]}_{res_num[1]}"
                        if all([altloc in self.ordered_atom_lookup.get_altlocs() for altloc in altlocs]): # XXX
                            clashes.append((name, res_num, badness, altlocs))
                return clashes

            unflipped_water_dict = {}
            for altloc in self.ordered_atom_lookup.get_altlocs():
                unflipped_water_dict[altloc]=altloc
            #clashes =  get_clashes(model_handle) 
            if constraint_weights[ConstraintsHandler.TwoAtomPenalty]>0:
                original_clashes = get_clashes(UntangleFunctions.model_handle(self.restrained_model_path))
                print(f"Original structure has {len(original_clashes)} clashes")
                if len(original_clashes) <= 5:
                    for clash in original_clashes:
                        print(clash)
            else:
                original_clashes=[]


        
        constraints_handler=ConstraintsHandler()
        constraints_to_skip=[]
        constraints_to_skip = [kind for kind,value in constraint_weights.items() if value <= 0]

        potential_overrides = None
        if ALLOW_OUTLIERS_FROM_POTENTIAL_OVERRIDES:
            potential_overrides = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}_potential_overrides.txt"
        constraints_handler.load_all_constraints(self.model_path,self.ordered_atom_lookup,symmetries=self.symmetries,water_water_nonbond=water_water_nonbond,constraints_to_skip=constraints_to_skip,
                                                 two_atom_penalty_tuples=original_clashes,
                                                 outliers_to_ignore_file=potential_overrides)
        print("Nonordered constraint properties loaded.")
        #for n,atom in enumerate(self.ordered_atom_lookup.select_atoms_by(names=["CA","C","N"])):

        print("Creating ordered atom chunks")

        
        atom_chunks:dict[str,AtomChunk] = {}
        def atom_id(atom:Atom):  # TODO Deprecated. Replace with OrderedTag
            def atom_id_from_params(atom_name,atom_altloc,atom_res_seq_num):
                return f"{atom_name}.{atom_altloc}{atom_res_seq_num}"
            #return f"{atom.get_name()}.{atom.get_altloc()}{OrderedAtomLookup.atom_res_seq_num(atom)}"
            return atom_id_from_params(atom.get_name(),atom.get_altloc(),OrderedAtomLookup.atom_res_seq_num(atom))
        
        ordered_atoms = self.ordered_atom_lookup.select_atoms_by()
        altloc_counts={}
        last_site=None
        for n,atom in enumerate(ordered_atoms):
            if atom.get_altloc() not in altloc_counts:
                altloc_counts[atom.get_altloc()]=0
            altloc_counts[atom.get_altloc()]+=1
            
            if n%1000==0:
                print(f"Creating chunk {n} / {len(ordered_atoms)} ")
            if atom.element=="H":
                continue
            is_water = UntangleFunctions.res_is_water(atom.get_parent())
            
            # Add alternate coord options for each disordered atom (coords and index of the option for each ordered atom)
            alt_pos_options:dict[int,Atom] = {}
            site_tag=DisorderedTag.from_atom(atom)
            if site_tag in self.ordered_atom_lookup.alt_pos_options:
                for option_index, alt_coords_dict in enumerate(self.ordered_atom_lookup.alt_pos_options[site_tag]):
                    if atom.get_altloc() in alt_coords_dict:
                        alt_pos_options[option_index+1]=alt_coords_dict[atom.get_altloc()]
                
            
            atom_chunks[atom_id(atom)]=AtomChunk(
                                         altloc_counts[atom.get_altloc()]+1,
                                         atom.get_altloc(),
                                         OrderedAtomLookup.atom_res_seq_num(atom),
                                         atom.get_parent(),
                                         atom,
                                         is_water,
                                         constraints_handler,
                                         alt_pos_options,
                                         )
            
            
        del ordered_atoms
        chunk_sets.append(atom_chunks)
        connection_types.append(LP_Input.Geomection)


        ############ Calculate wE for different combinations (don't need known constraints).
        '''
        possible_connections:dict[str,dict[str,LP_Input.ChunkConnection]]={}
        
        # Add chunk connections
        for chunk_set,connection_type in zip(chunk_sets,connection_types):
            for n,A in enumerate(chunk_set):
                possible_connections[A.unique_id()]={}
                for m,B in enumerate(chunk_set):
                    if B==A:
                        continue
                    if connection_type is LP_Input.Geomection:
                        if n>=m:
                            continue
                        connection = self.Geomection(A,B)
                        connection.calculate_distance()
                        if connection.ts_distance is not None:
                            #print(connection.ts_distance)
                            possible_connections[A.unique_id()][B.unique_id()]=connection
                    elif connection_type is LP_Input.ChunkConnection:
                        if B.start_resnum-A.end_resnum!=1:
                            continue
                        connection = self.ChunkConnection(A,B)
                        connection.calculate_distance(self.model_path,tmp_out_folder_path,self.ordered_atom_lookup.disordered_waters,quick_wE=quick_wE,dry=dry_run)
                        possible_connections[A.unique_id()][B.unique_id()]=connection
                        #print(connection.ts_distance)
                    else: assert False, connection_type
                for connection in possible_connections[A.unique_id()].values():
                    print(connection.A.get_disordered_tag(),connection.A.altloc,connection.B.get_disordered_tag(),connection.B.altloc,connection.ts_distance)
        '''
        ###########

        #################################################################
        print("Computing connection costs")
        possible_connections:list[LP_Input.Geomection]=[]
        # Apply when constraints don't involve water
        protein_constraints_that_include_H=[ConstraintsHandler.TwoAtomPenalty]  # Since the purpose of two atom penalty is to say "the current thing is wrong", in which case using the hydrogens is fine.
        # Apply when constraints involve water
        water_constraints_that_include_H=[ConstraintsHandler.TwoAtomPenalty,ConstraintsHandler.ClashConstraint,ConstraintsHandler.NonbondConstraint]
        for c, constraint in enumerate(constraints_handler.constraints):
            if c%1000 == 0: 
                print(f"Calculating constraint {c} / {len(constraints_handler.constraints)} ({constraint}) ")

            if any([tag.resnum() in self.ordered_atom_lookup.water_residue_nums for tag in constraint.site_tags]):
                constraints_that_include_H=water_constraints_that_include_H
            else:
                constraints_that_include_H=protein_constraints_that_include_H
            atoms_for_constraint = self.ordered_atom_lookup.select_atoms_for_sites(
                constraint.site_tags,
                exclude_H=type(constraint) not in constraints_that_include_H 
            )
            res_name_num_dict={}
            for a in atoms_for_constraint:
                res_name_num_dict[a] = (a.name,OrderedAtomLookup.atom_res_seq_num(a) ) 
            # Generate each combination of n atoms
            combinations_iterator:itertools.combinations[tuple[Atom]] = itertools.combinations(atoms_for_constraint, constraint.num_atoms())
            
            debug_print=False
            if debug_print:
                combinations_iterator = list(combinations_iterator)
                print(len(list(combinations_iterator)))
                old_num_connections=len(list(possible_connections))
            # print(set([res_name_num_dict[a] for a in atoms]))
            # print(atoms)
            # assert False
            # NOTE that the tension is contributed by ALL geometric restraints affecting the sites, not just the one represented by 'constraint'
            tension_mod=None
            if self.APPLY_TENSION_MOD is not None:
                constraint.set_tension(np.sum([self.APPLY_TENSION_MOD[tag] for tag in constraint.site_tags])/len(constraint.site_tags))
                constraint.set_max_site_tension(max([self.APPLY_TENSION_MOD[tag] for tag in constraint.site_tags]))
                if APPLY_TENSION_MOD:
                    #tension_mod = (constraint.get_tension()+1)**2
                    #tension_mod = (constraint.get_max_site_tension()+1)**2
                    tension_mod = (constraint.get_max_site_tension()+1)**2
                    assert tension_mod>=1


            for atoms in combinations_iterator:
                # Don't have two alt locs of same atom in group
                if len(set([res_name_num_dict[a] for a in atoms]))!= len(atoms):                     
                    continue
                output = constraint.get_cost(atoms,scoring_function)
                if output is not None:
                    ideal,z_score,distance = output
                    weight_mod=1
                    weight_mod*=constraint_weights[type(constraint)]

                    # print(constraint.atom_id)
                    # print(constraint.atom_names(),constraint.residues())
                    
                    atom_chunks_selection:list[AtomChunk] = []
                    impossible_H_parent=False
                    contains_H=False
                    for a in atoms:
                        if a.get_name()[0] != "H":
                            atom_chunks_selection.append(atom_chunks[atom_id(a)])
                        else:
                            contains_H=True
                            # Convert restraints on riding hydrogens to restraints on parent atoms
                            assert type(constraint) in constraints_that_include_H, (type(constraint),[a.get_name() for a in atoms])
                            res_num, altloc = OrderedAtomLookup.atom_res_seq_num(a), a.get_altloc()
                            
                            res_atoms = []
                            for name in self.ordered_atom_lookup.better_dict[res_num]:
                                res_atoms.extend(self.ordered_atom_lookup.better_dict[res_num][name].values())
                            possible_parents: list[Atom] = self.ordered_atom_lookup.select_atoms_from(
                                res_atoms,
                                altlocs=[altloc],
                                exclude_H=True,
                            )
                            parent_name = UntangleFunctions.H_get_parent_fullname(
                                a.get_name(),
                                [p.get_fullname() for p in possible_parents]
                            ).strip()
                            parent_atom = self.ordered_atom_lookup.better_dict[res_num][parent_name][altloc]
                            # Commenting out for now because breaks plotting code. #TODO should work now but need to test get same outcome.
                            # if parent_atom in atoms:
                            #     for p_a in self.ordered_atom_lookup.better_dict[res_num][parent_name]:
                            #         if p_a != parent_atom:
                            #             impossible_H_parent=True
                            #             break
                            #     if impossible_H_parent:
                            #         break


                            #assert len(parent_atom) == 1, (parent_atom,parent_name,res_num,altloc)
                            #parent_atom = parent_atom[0]
                            atom_chunks[atom_id(parent_atom)].constraints_holder.add(constraint,None) # NOTE XXX
                            atom_chunks_selection.append(atom_chunks[atom_id(parent_atom)])
                    


                    
                    
                    
                    # if impossible_H_parent:
                    #     continue # Commenting out for now because breaks plotting code

                    num_hydrogens = len([a.get_name() for a in atoms if a.element=="H"])
                    hydrogens=None
                    if num_hydrogens > 0:
                        hydrogens:list[str] = [a.get_name() if a.element=="H" else "x" for a in atoms]
                    
                    for ach in atom_chunks_selection:
                        assert constraint in ach.constraints_holder.constraints, (constraint, ach.get_disordered_tag(), [c for c in ach.constraints_holder.constraints])

                    if tension_mod is not None:
                        weight_mod*=tension_mod
                    if constraint.outlier_ok:
                        weight_mod*=weight_mod_for_allowed_outliers
                    if contains_H:
                        weight_mod*=hydrogen_restraint_scale



                    z_scores=[z_score]
                    distances=[distance]
                    position_option_indices_list=[(0,)*len(atom_chunks_selection)]

                    # Duct-tape in alternate costs/constraints/geomections for alternate positions XXX
                    # XXX redundancy, atom chunks have access to atoms.
                    
                    have_alternate_positions=False
                    for ch in atom_chunks_selection:
                        if len(ch.alt_pos_options)>0:
                            have_alternate_positions=True
                            break
                    if have_alternate_positions:
                        idx_altatom_tuples = [[(0,a)]+ [(i,alt_atoms) for i,alt_atoms in ch.alt_pos_options.items()] for a, ch in zip(atoms,atom_chunks_selection)]
                        alternate_combinations = itertools.product(*idx_altatom_tuples)
                        for alt_comb in alternate_combinations:
                            alt_position_indices, alt_atoms = zip(*alt_comb)
                            if all([i==0 for i in alt_position_indices]):
                                # Non-alternate position - already calculated. # XXX 
                                continue
                            assert len(alt_position_indices)==len(position_option_indices_list[0])
                            position_option_indices_list.append(alt_position_indices)
                            alt_ideal,alt_z,alt_d = constraint.get_cost(alt_atoms,scoring_function)
                            assert alt_ideal==ideal
                            distances.append(alt_d)
                            z_scores.append(alt_z)
                    
                    for d,z,pos_indices in zip(distances,z_scores,position_option_indices_list):
                        if any([i != 0 for i in pos_indices]):
                            d*=ALTERNATE_POSITIONS_FACTOR
                            z*=ALTERNATE_POSITIONS_Z_SCORE_FACTOR
                        connection = self.Geomection(atom_chunks_selection,d*weight_mod,pos_indices,
                                                     type(constraint),hydrogens,ideal,z,
                                                     constraint.get_max_site_tension(),constraint.outlier_ok)  # XXX putting max site tension in here is bad
                        possible_connections.append(connection)
            if debug_print:
                print("added:",len(list(possible_connections))-old_num_connections)
                print("total:",len(list(possible_connections)))
                        #print([(ch.name,ch.resnum,distance) for ch in connection.atom_chunks])
        #################################################################

        

        disordered_connections:dict[str,list[LP_Input.Geomection]] ={} # options for each alt connection
        
        
        for connection in possible_connections:
            connection_id = connection.get_disordered_connection_id()
            if connection_id not in disordered_connections:
                disordered_connections[connection_id]=[]
            assert connection not in disordered_connections[connection_id] 
            for other in disordered_connections[connection_id]:
                assert (connection.atom_chunks!=other.atom_chunks) or (connection.hydrogen_tag!=other.hydrogen_tag) \
                or (connection.connection_type!=other.connection_type) or (connection.position_option_indices != other.position_option_indices) , f"duplicate connections of kind {connection.connection_type}, involving ordered atoms {[ch.unique_id() for ch in connection.atom_chunks]}, with badness {connection.ts_distance} and {other.ts_distance}!"
            disordered_connections[connection_id].append(connection)
    
        if weight_for_range:
            for disordered_connection_id, ordered_connections in disordered_connections.items():
                z_scores = [conn.z_score for conn in ordered_connections if conn.single_altloc()]
                min_z = min(z_scores)                        
                max_z = max(z_scores)   

                range_threshold=0.7
                max_z_threshold=1.3
                if (max_z - min_z)>range_threshold and max_z > max_z_threshold:
                    for conn in ordered_connections:
                        conn.ts_distance*=10*(1+(max_z-min_z))**2
        
        if improvement_factor !=1:
            assert False, "Not implemented"
            for disordered_connection_id, ordered_connections in disordered_connections.items():
                for conn in disordered_connection_id:
                    pass
        if suppress_worse:
            for disordered_connection_id, ordered_connections in disordered_connections.items():
                max_dist  = max([conn.ts_distance for conn in ordered_connections if conn.single_altloc()])

                for conn in ordered_connections:
                    if conn.ts_distance > max_dist:
                        conn.ts_distance = max_dist + (conn.ts_distance-max_dist)/2 #np.log(1+conn.ts_distance-max_dist) 
                        
        if force_solution_reference is not None:
            for disordered_connection_id, ordered_connections in disordered_connections.items():
                for conn in ordered_connections:
                    ordered_tags = [ach.get_ordered_tag() for ach in conn.atom_chunks]
                    solution_to_altlocs = [force_solution_reference[tag] for tag in ordered_tags]
                    # If connection is between two atoms that are not in the same conformation in the forced solution, forbid it. 
                    conn.forbidden = len(set(solution_to_altlocs))!=1
                    # if ordered_tags[0].resnum() <= 8 and len(ordered_tags)==2:
                    #     print(ordered_tags, "ON" if not conn.forbidden else "")
                    #     print([ach.altloc for ach in conn.atom_chunks],solution_to_altlocs)
                    
                # forbid connections between
        else:
            # If the current connections have a tension above a certain value, let solver consider all alternatives.
            print(f"Re-enabling connections that are alternatives to current connections with tension > {LP_Input.min_tension_where_anything_goes}")
            num_bad_current_disordered_connections={k:0 for k in LP_Input.min_tension_where_anything_goes}
            num_connections_re_enabled={k:0 for k in LP_Input.min_tension_where_anything_goes}
            for disordered_connection_id, ordered_connections in disordered_connections.items():

                conn_type = ordered_connections[0].connection_type
                if conn_type not in LP_Input.min_tension_where_anything_goes:
                    continue
                if ordered_connections[0].max_site_tension >= LP_Input.min_tension_where_anything_goes[conn_type]:
                    num_bad_current_disordered_connections[conn_type]+=1
                else: continue

                for conn in ordered_connections:
                    if conn.forbidden:
                        num_connections_re_enabled[conn.connection_type]+=1
                    conn.forbidden=False
                    if conn.single_altloc():
                        if conn_type == ConstraintsHandler.AngleConstraint:
                            conn.ts_distance*=high_tension_penalty
            print(f"Number of high tension disordered connections detected: {num_bad_current_disordered_connections}") 
            print(f"Ordered connections re-enabled: {num_connections_re_enabled}")
            # If the current connections have a sigma above a certain value, let solver consider all alternatives.

            if not TURN_OFF_MIN_SIGMAS:
                print(f"Re-enabling connections that are alternatives to current connections with sigma > {LP_Input.min_sigmas_where_anything_goes}")
                num_bad_current_disordered_connections={k:0 for k in LP_Input.min_sigmas_where_anything_goes}
                num_connections_re_enabled={k:0 for k in LP_Input.min_sigmas_where_anything_goes}
                for disordered_connection_id, ordered_connections in disordered_connections.items():
                    if ordered_connections[0].outlier_ok:  # XXX represents disordered connection
                        continue
                    
                    if ordered_connections[0].connection_type not in LP_Input.min_sigmas_where_anything_goes:
                        continue
                    for conn in ordered_connections:
                        if conn.single_altloc() and conn.z_score >= LP_Input.min_sigmas_where_anything_goes[conn.connection_type]:
                            num_bad_current_disordered_connections[conn.connection_type]+=1
                            break
                    else: continue
                    for conn in ordered_connections:
                        if conn.forbidden:
                            num_connections_re_enabled[conn.connection_type]+=1
                        conn.forbidden=False
                print(f"Number of high sigma disordered connections detected: {num_bad_current_disordered_connections}") 
                print(f"Ordered connections re-enabled: {num_connections_re_enabled}")







                


        #finest_depth_chunks=orderedResidues
        finest_depth_chunks=list(atom_chunks.values())
        return finest_depth_chunks,disordered_connections

    @staticmethod
    def create_altloc_subset_model(model,altloc_subset):
        struct=PDBParser().get_structure("struct",model)
        atom_lookup = OrderedAtomLookup(struct.get_atoms(),
                                                    protein=True,waters=True,
                                                    altloc_subset=altloc_subset)   
        #temp_path=working_model[:-4]+"_subsetOut.pdb"
        temp_path=LP_Input.subset_model_path(model,altloc_subset)[:-4]+"Out.pdb"
        atom_lookup.output_as_pdb_file(reference_pdb_file=model,out_path=temp_path)
        return temp_path