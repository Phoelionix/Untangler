#%%
from typing import Any
from LinearOptimizer import Solver
from LinearOptimizer.Input import OrderedAtomLookup, ConstraintsHandler
from LinearOptimizer.Swapper import Swapper
import UntangleFunctions
from UntangleFunctions import assess_geometry_wE, get_R,pdb_data_dir,create_score_file,get_score,score_file_name,geo_file_name,res_is_water, parse_symmetries_from_pdb
import subprocess
import os, sys
import numpy as np
import shutil
import random
from types import SimpleNamespace 
from multiprocessing import Pool
from time import sleep
from enum import Enum
import itertools
from Bio.PDB import PDBParser,Structure,PDBIO
from Bio.PDB.Atom import Atom,DisorderedAtom
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KNeighborsClassifier
import psutil
from LinearOptimizer.Tension import GeoXrayTension
from matplotlib import pyplot as plt
#import create_delta_constraints_file



DISABLE_WATER_ALTLOC_OPTIM=True
TURN_OFF_BULK_SOLVENT=False
CONSIDER_WE_WHEN_CHOOSING_BEST_BATCH=True
PHENIX_ORDERED_SOLVENT=False
TENSIONS=False  # Enables behaviours of re-enabling connection options involving high-tension sites, and, if option enabled, to scale cost by tensions
PHENIX_FREEZE_WATER=False
PHENIX_DISABLE_CDL=False # TODO -> set True
PROPOSE_IGNORES_H=False

assert not PROPOSE_IGNORES_H, "Comparison with previous best not implemented properly"

# TODO:
# down weight swaps that were previously made but need to be made again (i.e. cases where it's not tangled and the density is pushing it a different way.)
# Try forbid connection changes in sidechain. 
# Make flag in swap options that says what the altlocs were swapped around. So that on debug reruns the same altlocs can be used.

class Untangler():
    working_dir = os.path.abspath(os.getcwd())
    output_dir = f"{working_dir}/output/"
    refine_shell_file=f"Refinement/Refine.sh"
    refine_refmac_shell_file=f"Refinement/Refine_refmac.sh"
    ####
    # Skip stages. Requires files to already have been generated (up to the point you are debugging).
    # If a model file hasn't been generated, will use the latest filename that exists and is expected to have been generated up to that point in the process.
    debug_skip_refine = False  # Note: Can set to True alongside debug_skip_first_swaps to skip to first proposal
    debug_phenix_ordered_solvent_on_initial=False
    debug_skip_initial_refine=True
    debug_skip_first_unrestrained_refine=False
    never_do_unrestrained=False
    debug_skip_first_swaps=False
    debug_skip_first_focus_swaps=False # many swaps strategy only 
    debug_skip_unrestrained_refine=False
    debug_skip_holton_data_generation=False
    debug_always_accept_proposed_model=False
    auto_group_waters=False
    debug_skip_to_loop=0
    debug_skip_initial_holton_data_generation =debug_skip_initial_refine or (debug_skip_to_loop!=0)
    refmac_refine_water_occupancies_initial=False
    ##
    PHENIX = 1
    REFMAC = 2
    refinement=PHENIX
    ##
    O_bond_change_period=5
    ####
    num_threads=15
    class Score():
        def __init__(self,combined,wE,R_work,R_free):
            self.wE=wE
            self.combined=combined  
            self.R_work=R_work
            self.R_free = R_free
        def __repr__(self):
            return f"{self.combined} | wE: {self.wE} Rwork: {self.R_work} Rfree: {self.R_free}"
        @staticmethod
        def inf_bad_score()->'Untangler.Score':
            return Untangler.Score(np.inf,np.inf,np.inf,np.inf)

    # TODO keep refining while Rfree decreasing.
    def __init__(self,acceptance_temperature=1,max_wE_frac_increase=0, num_end_loop_refine_cycles=8,  #8,
                 default_wc=1,wc_anneal_start=1,wc_anneal_loops=0, starting_num_best_swaps_considered=20, # 20,
                 max_num_best_swaps_considered=100,num_refine_for_positions_macro_cycles_phenix=3,
                 num_loops_water_held=0,weight_factors=None,
                 max_bond_changes=9999,altloc_subset_size=3,refine_for_positions_geo_weight=0.1):
        self.set_hyper_params(acceptance_temperature,max_wE_frac_increase,num_end_loop_refine_cycles,
                              default_wc,wc_anneal_start,wc_anneal_loops, starting_num_best_swaps_considered,
                              max_num_best_swaps_considered,num_refine_for_positions_macro_cycles_phenix, 
                              num_loops_water_held,max_bond_changes,altloc_subset_size,refine_for_positions_geo_weight)
        self.previously_swapped = []
        self.model_handle=None
        self.current_model=None
        self.initial_score=self.current_score=self.best_score=None
        self.best_score=self.current_score
        self.swapper:Swapper=None
        self.swaps_history:list[Swapper.SwapGroup] = [] 
        self.model_protein_altlocs=None
        self.model_solvent_altlocs=None
        self.weight_factors = weight_factors
        
        os.makedirs(UntangleFunctions.separated_conformer_pdb_dir(),exist_ok=True)
    def set_hyper_params(self,acceptance_temperature=1,max_wE_frac_increase=0, num_end_loop_refine_cycles=2,  # 8,
                 default_wc=1,wc_anneal_start=1,wc_anneal_loops=0, starting_num_best_swaps_considered=5, # 20,
                 max_num_best_swaps_considered=100,num_refine_for_positions_macro_cycles_phenix=3,
                 num_loops_water_held=0,
                 max_bond_changes=None,altloc_subset_size=3,refine_for_positions_geo_weight=0.1):
        # NOTE Currently max wE increase is percentage based, 
        # but TODO optimal method needs to be investigated.
        self.n_best_swap_start=starting_num_best_swaps_considered
        self.n_best_swap_max = max_num_best_swaps_considered
        self.num_refine_for_positions_macro_cycles_phenix=num_refine_for_positions_macro_cycles_phenix
        self.max_wE_frac_increase=max_wE_frac_increase
        self.acceptance_temperature=acceptance_temperature
        self.num_end_loop_refine_cycles=num_end_loop_refine_cycles
        self.num_loops_water_held=num_loops_water_held
        self.default_wc=default_wc
        self.wc_anneal_start = wc_anneal_start
        self.wc_anneal_loops=wc_anneal_loops
        self.max_bond_changes=max_bond_changes
        self.altloc_subset_size=altloc_subset_size
        self.refine_for_positions_geo_weight=refine_for_positions_geo_weight

    def delete_zero_occupancy_waters(self,pdb_path,out_path):
        assert os.path.abspath(pdb_path) != os.path.abspath(out_path)
        solvent_res_names=["HOH"]
        start_strs_considered = ["ATOM","HETATM"]
        def replace_res_num(line,res_num):
            res_num = str(res_num)
            res_num = ' '*(4-len(res_num))+res_num
            return line[:22]+res_num+line[26:]
        def replace_serial_num(line,serial_num):
            serial_num = str(serial_num)
            serial_num = ' '*(5-len(serial_num))+serial_num
            return line[:6]+serial_num+line[11:]
        with open(pdb_path) as I, open(out_path,'w') as O:
            max_resnum=0
            max_serial_num=0
            for line in I:
                for s in start_strs_considered:
                    if line.startswith(s):
                        resname = line[17:20]
                        if resname not in solvent_res_names:
                            resnum = line[22:26]
                            serial_num = line[6:11]
                            max_resnum=max(max_resnum,int(resnum))
                            max_serial_num=max(max_serial_num,int(serial_num))
                        else:
                            break
                else: # Not modifying
                    O.write(line)
                    continue
                occ = float(line[54:60].strip())
                if occ == 0:
                    continue
                else:
                    max_resnum+=1
                    max_serial_num+=1
                    line = replace_res_num(line,max_resnum)
                    line = replace_serial_num(line,max_serial_num)
                    O.write(line)
                    

    def prepare_pdb_and_read_altlocs(self,pdb_path,out_path,sep_chain_format=False,altloc_from_chain_fix=False):
        # Gets into format we expect. !!!!!!Assumes single chain!!!!!
        def replace_occupancy(line,occ):
            occ=f"{occ:.3f}"
            occ = ' '*(6-len(occ))+occ
            return line[:54] + occ + line[60:]
        def replace_chain(line,chain_id):
            chain_id = str(chain_id)
            assert len(chain_id)==1
            return line[:21]+chain_id+line[22:]
        def replace_serial_num(line,serial_num):
            serial_num = str(serial_num)
            serial_num = ' '*(5-len(serial_num))+serial_num
            return line[:6]+serial_num+line[11:]
        def replace_altloc(line,altloc):
            return line[:16]+altloc+line[17:]
        def replace_res_num(line,res_num):
            res_num = str(res_num)
            res_num = ' '*(4-len(res_num))+res_num
            return line[:22]+res_num+line[26:]
            
        protein_altlocs = []
        solvent_altlocs = []
        with open(pdb_path) as I:
            max_resnum=0
            start_lines = []
            solvent_lines=[]
            end_lines = []
            atom_dict:dict[str,dict[str,dict[str,str]]] = {}  
            last_chain=None
            solvent_res_names=["HOH"]
            solvent_chain_id = "z"
            warned_collapse=False

            for line in I:
                if line.startswith("TER") or line.startswith("ANISOU"):
                    continue
                start_strs_considered = ["ATOM","HETATM"]
                for s in start_strs_considered:
                    if line.startswith(s):
                        break
                else: # Not modifying
                    if len(atom_dict)==0:
                        start_lines+=line
                    else:
                        end_lines += line
                    continue

                # Modifying
                name = line[12:16].strip()
                altloc = line[16]
                resname = line[17:20]
                space = line[20]
                chain = line[21]
                resnum = int(line[22:26])
                if altloc == ' ' and altloc_from_chain_fix:
                    altloc = chain
                    line = replace_altloc(line,altloc)
                if resname in solvent_res_names:
                    solvent_lines.append(replace_chain(line,solvent_chain_id))
                    if altloc not in solvent_altlocs:
                        solvent_altlocs.append(altloc) 
                    continue
                assert len(end_lines)==0
                
                # Non-solvent atoms
                if not sep_chain_format and not warned_collapse and chain != last_chain and last_chain is not None:
                    print("Warning: Multiple chains detected. Collapsing chains into single chain")
                    warned_collapse=True
                if resnum not in atom_dict:
                    atom_dict[resnum] = {}
                
                assert (altloc != ' '), line 

                if altloc not in atom_dict[resnum]:
                    atom_dict[resnum][altloc] = {}
                    
                    if altloc not in protein_altlocs:
                        protein_altlocs.append(altloc) 
                
                atom_dict[resnum][altloc][name]=line  
                max_resnum=max(resnum,max_resnum)
                last_chain = chain
                continue
                    
        n=0
        # Add non-solvent atoms
        if not sep_chain_format: # format for untangler stuff
            protein_chain_id = "A"
            for res_atom_dict in atom_dict.values():
                for altloc_atom_dict in res_atom_dict.values():
                    for line in altloc_atom_dict.values():
                        n+=1
                        modified_line = line
                        
                        modified_line = replace_occupancy(modified_line,
                            1/len(protein_altlocs)) # Set occupancies to all be same
                        modified_line = replace_chain(modified_line,protein_chain_id)
                        modified_line = replace_serial_num(modified_line,n)
                        start_lines.append(modified_line)
        else: # Note that lines for each chain need to be contiguous in the file
            chain_dict={}
            for res_atom_dict in atom_dict.values():
                for altloc, altloc_atom_dict in res_atom_dict.items():
                    protein_chain_id = altloc
                    for line in altloc_atom_dict.values():
                        n+=1
                        modified_line = line
                        modified_line = replace_occupancy(modified_line,
                            1/len(protein_altlocs)) # Set occupancies to all be same
                        modified_line = replace_chain(modified_line,protein_chain_id)
                        modified_line = replace_serial_num(modified_line,n)
                        if altloc not in chain_dict:
                            chain_dict[altloc]=[]
                        chain_dict[altloc].append(modified_line)
            for _, lines in chain_dict.items():
                for modified_line in lines:
                    start_lines.append(modified_line)


        # # Make sure waters don't share residue numbers with protein
        # for line in solvent_lines:
        #     max_resnum+=1
        #     modified_line = replace_res_num(line,max_resnum)
        #     start_lines.append(modified_line)

        # Make sure waters don't share residue numbers with protein
        min_solvent_resnum=99999999
        for line in solvent_lines:
            solvent_resnum=int(line[22:26])
            min_solvent_resnum = min(solvent_resnum,min_solvent_resnum)
        shift = max_resnum-min_solvent_resnum + 1
        for line in solvent_lines:
            solvent_resnum=int(line[22:26])
            # in case of gaps...
            assert (solvent_resnum+shift) - max_resnum >=0, (max_resnum,solvent_resnum,shift)
            if (solvent_resnum+shift)-max_resnum > 1:
                shift = max_resnum-solvent_resnum + 1 
            max_resnum = shift+solvent_resnum
            
            modified_line = replace_res_num(line,max_resnum)
            start_lines.append(modified_line)

        with open(out_path,'w') as O:
            O.writelines(start_lines+end_lines)
        self.model_protein_altlocs=set(protein_altlocs)
        self.model_solvent_altlocs=set(solvent_altlocs)

    def run(self,pdb_file_path,hkl_file_path,desired_score=18.6,max_num_runs=100):
        # pdb_file_path: path to starting model
        # hkl_file_path: path to reflection data.
        # TODO Currently assume in data folder.
        assert hkl_file_path[-4:]==".mtz", f"hkl path doesn't end in '.mtz': {hkl_file_path}"
        self.hkl_path = os.path.abspath(hkl_file_path) 
        #assert os.path.dirname(hkl_file_path)[-5:-1]=="data", hkl_file_path    
        assert pdb_file_path[-4:]==".pdb", f"file path doesn't end in '.pdb': {pdb_file_path}"
        self.model_handle = os.path.basename(pdb_file_path)[:-4]
        #self.current_model= Untangler.output_dir+f"{self.model_handle}_current.pdb"
        self.current_model = Untangler.output_dir+f"{self.model_handle}_start.pdb"
        os.makedirs(Untangler.output_dir,exist_ok=True)
        output_archive = Untangler.output_dir+".archive/"
        os.makedirs(output_archive,exist_ok=True)

        move_all_to_archive=False
        for p in os.listdir(Untangler.output_dir):
            fixed_folders=["refine_logs"]  # TODO put optimizer logs xLO in here
            fp = Untangler.output_dir+p
            if os.path.abspath(fp) != os.path.abspath(output_archive):
                if p in fixed_folders:
                    for sub_p in os.listdir(fp):
                        sub_file_path=os.path.join(fp,sub_p)
                        assert os.path.isfile(sub_file_path)  # Okay to turn this off, just at the moment don't expect fixed folders (refine_logs) to have subfolders.
                        os.makedirs(os.path.join(output_archive,p,""),exist_ok=True)
                        shutil.move(sub_file_path,os.path.join(output_archive,p,sub_p))
                else:    
                    if move_all_to_archive:
                        shutil.move(fp,output_archive+p)

        #shutil.copy(pdb_file_path,self.current_model)
        self.prepare_pdb_and_read_altlocs(pdb_file_path,self.current_model)
        self.symmetries = parse_symmetries_from_pdb(pdb_file_path) #TODO do in LinearOptimizer/Input. But first need to figure out why REMARK 290 is being removed... or work around it.
        
        if self.auto_group_waters:
            if self.model_solvent_altlocs != self.model_protein_altlocs:
                self.group_waters_to_altlocs(self.current_model,self.current_model)
            assert self.model_solvent_altlocs==self.model_protein_altlocs
        self.swapper = Swapper()


        self.first_loop=self.loop=0 if self.debug_skip_to_loop is None else self.debug_skip_to_loop
        
        skip_init = (self.debug_skip_refine or self.debug_skip_initial_refine or self.first_loop>0)
        if not skip_init:
            initial_model=self.initial_refine(self.current_model,debug_skip=skip_init) 
        else:
            initial_model=self.current_model
            if self.first_loop>0:
                initial_model= self.current_model_name_at_loop(self.first_loop)
        assert os.path.exists(initial_model), f"expected model file not found at {initial_model}"
        
        self.current_model=initial_model

        if not self.debug_skip_holton_data_generation and not self.debug_skip_initial_holton_data_generation:
            self.initial_score = Untangler.Score(*assess_geometry_wE(initial_model, log_out_folder_path=self.output_dir))
        else:
            if not os.path.exists(score_file_name(initial_model)):
                create_score_file(initial_model,log_out_folder_path=self.output_dir)
            self.initial_score = Untangler.Score(*get_score(score_file_name(initial_model)))
        self.current_score = self.best_score = Untangler.Score.inf_bad_score()# i.e. force accept first solution
        if self.first_loop>0:
            self.current_score=self.best_score=self.initial_score
        # if initial_model!=self.current_model:
        #     shutil.copy(initial_model,self.current_model)
        
        while (self.current_score.combined > desired_score and self.loop < max_num_runs): 
            print(f"Initial score: {self.initial_score}")
            print(f"Best score: {self.best_score}")
            print(f"Current score: {self.current_score}")
            self.step()
            self.loop+=1

        if self.current_score.combined > desired_score:
            print(f"Failed to reach target score of {desired_score}")
        else:
            print("Target score reached")
            print("Initial:", self.initial_score)
            print("Best:",self.best_score)
            print("Final:",self.current_score)

        print(self.current_score)


    def step(self):
        self.refinement_loop()


    def skip_swaps(self):
        return self.debug_skip_first_swaps and self.loop == self.first_loop
    def swap_specific_residues(self,resnames,swapper,tensions,model_to_swap:str,file_tag,altloc_subset_size=15,cycles=1):
        if not self.skip_swaps():
            print(f"Optimizing connections of {', '.join(resnames)} residues")
        return self.many_swapped(swapper,tensions,model_to_swap,True,altloc_subset_size=altloc_subset_size,allowed_resnames=resnames,cycles=cycles,file_tag=file_tag,
                                 forbidden_atom_bond_changes=["C","N"],  # Because otherwise could ruin angle with not considered atoms. TODO need to automate this sort of decision-making! If only some atoms are allowed, need to include all restraints that involve the atoms even those that do not involve allowed atoms
                                 forbid_altloc_changes=["C","N"]) # Bad (as in highly likely missing better solutions), but currently necessary for solve time to be feasible. Presumably because haven't included connections that involve "non-active" atoms
    def swap_cysteines(self,swapper,tensions,model_to_swap:str,altloc_subset_size=15):
       return self.swap_specific_residues(["CYS"],swapper,tensions,model_to_swap,file_tag="CysSwaps",altloc_subset_size=altloc_subset_size)


    def get_altloc_subsets(self,altloc_subset_size,num_combinations):
        altloc_subsets = [None]
        if len(self.model_protein_altlocs)>altloc_subset_size:
            altloc_subsets = []
            altlocs_tmp = [altloc for altloc in self.model_protein_altlocs]
            random.shuffle(altlocs_tmp)
            assert ' ' not in self.model_protein_altlocs
            while len(altlocs_tmp) >= altloc_subset_size and len(altloc_subsets) < num_combinations:  
                altloc_subsets.append(altlocs_tmp[:altloc_subset_size])
                del altlocs_tmp[:altloc_subset_size]
        return altloc_subsets

    #TODO resample if sample same subsets as last cycle
    def many_swapped(self,swapper,tensions,model_to_swap:str,restrained_refine_pdb_file_path:str,allot_protein_independent_of_waters:bool,altloc_subset_size=3,num_combinations=30,cycles=3,conformer_stats=False,
                     file_tag="manySwaps", allowed_resnums=None,allowed_resnames=None,forbidden_atom_bond_changes=[],forbidden_atom_any_connection_changes=[],forbid_altloc_changes=[]):
        
        # forbidden_atom_bond_changes:  atom names for which  bonds involving them should not be changed.
        
        #TODO try strategy of making one altloc as good as possible, while other can be terrible.
        
        #TODO try strategy of focusing on worst conformers

        working_model = f"{self.output_dir}/{self.model_handle}_{file_tag}{self.loop}.pdb"
        
        all_swaps=[]
        if self.skip_swaps():
            return working_model, ["Unknown"]

        
        measure_wE_after=False
        debug_prev_subset = None
        #debug_prev_subset=["D","B"]
        #debug_prev_subset=[ "A", "B"]
        #debug_prev_subset=[ "E", "C"]
        #debug_prev_subset=[ "D", "H"]
        #debug_prev_subset=[ "V", "B"]
        #debug_prev_subset=[ "t", "e","O"]
        #debug_prev_subset=["p", "C", "t", "J"]
        if debug_prev_subset is None or (not os.path.exists(working_model)):
            self.delete_zero_occupancy_waters(model_to_swap,working_model) 
            #shutil.copy(model_to_swap,working_model)
        if len(self.model_protein_altlocs)<=altloc_subset_size:
            cycles = 1 # No point repeating if considering all at once.


        if conformer_stats:
            Solver.LP_Input.prepare_geom_files(model_to_swap,self.model_protein_altlocs,water_swaps=False)

        for r in range(cycles):
            #altloc_subset_combinations:itertools.combinations[tuple[str]] = list(itertools.combinations(self.model_protein_altlocs, altloc_subset_size))
            #print("num possible combinations:",len(altloc_subset_combinations))
            #if len(altloc_subset_combinations) > num_combinations:
                #altloc_subset_combinations = random.sample(altloc_subset_combinations,num_combinations)

            # ITerate over unique sets, since we are creating geom files for all altloc subsets in the below loop
            altloc_subsets=self.get_altloc_subsets(altloc_subset_size,num_combinations)

            debug_skip_geom_file_prep=False
            #altloc_subsets=["teO"]
            #altloc_subsets=["tsO"]
            if debug_prev_subset is not None:
                measure_wE_after=True
                altloc_subsets=[debug_prev_subset]
            elif not debug_skip_geom_file_prep:
                UntangleFunctions.clear_geo() # Mainly space concerns.
                # TODO should be preparing geo file for one conformation. Barring clashes (and maybe some nonbonds?) it is all calculated in LinearOptimizer.input 
                # TODO swaps can create nonbond issues that are not recorded due to not being present in geo file?
                Solver.LP_Input.prepare_geom_files(working_model,altloc_subsets,allowed_resnames=allowed_resnames,
                                                      water_swaps=(altloc_subset_size==2))
            else:
                print("Warning: reusing old geom files")

            for n, altloc_subset in enumerate(altloc_subsets): 
                header=f"=========Altloc Allotment {n+1}/{len(altloc_subsets)}, Cycle {r+1}/{cycles}=========="
                print(f"\n{header}\n{'-'*len(header)}")
                if altloc_subset is None:
                    print("Optimizing connections")
                else:
                    print(f"Optimizing connections across altlocs {', '.join(altloc_subset)}")
                process = psutil.Process()
                print("Memory usage (MB):",psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)

                cand_models, cand_swaps = self.candidate_models_from_swapper(swapper,1,tensions,working_model,restrained_refine_pdb_file_path,allot_protein_independent_of_waters,
                                                                   altloc_subset=altloc_subset,need_to_prepare_geom_files=False,
                                                                   resnums=allowed_resnums,resnames=allowed_resnames,forbidden_atom_bond_changes=forbidden_atom_bond_changes,
                                                                   forbidden_atom_any_connection_changes=forbidden_atom_any_connection_changes,
                                                                   forbid_altloc_changes=forbid_altloc_changes)
                assert len(cand_models)==len(cand_swaps)==1
                if not debug_prev_subset:
                    shutil.move(cand_models[0],working_model)
                else:
                    shutil.move(cand_models[0],working_model[:-4]+"_debug_"+''.join(altloc_subset)+".pdb")
                all_swaps.extend(cand_swaps[0])

                if measure_wE_after:
                    struct=PDBParser().get_structure("struct",working_model)
                    ordered_atom_lookup = OrderedAtomLookup(struct.get_atoms(),
                                                                protein=True,waters=not allot_protein_independent_of_waters,
                                                                altloc_subset=altloc_subset)   
                    #temp_path=working_model[:-4]+"_subsetOut.pdb"
                    temp_path=Solver.LP_Input.subset_model_path(working_model,altloc_subset)[:-4]+"Out.pdb"
                    ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=working_model,out_path=temp_path)
                    assess_geometry_wE(temp_path,Solver.LP_Input.geo_log_out_folder())

                    if conformer_stats:
                        Solver.LP_Input.prepare_geom_files(working_model,self.model_protein_altlocs,water_swaps=False)

            if debug_prev_subset is not None:
                raise Exception("End debug")

        print("=======End Altloc Allotment=============\n")
        return working_model, all_swaps
    def candidate_models_from_swapper(self,swapper:Swapper,num_solutions,tensions,model_to_swap:str,restrained_refine_pdb_file_path:str,allot_protein_independent_of_waters:bool,altloc_subset=None,need_to_prepare_geom_files=True,read_prior_run=False,
                                      resnums=None,resnames=None,forbidden_atom_bond_changes=[],forbidden_atom_any_connection_changes=[],forbid_altloc_changes=[],
                                      MAIN_CHAIN_ONLY=False,SIDE_CHAIN_ONLY=False,NO_CB_CHANGES=False,num_sols_already_saved_this_loop=0): #TODO refactor as method of Swapper class
        # TODO should be running solver for altloc set partitioned into subsets, not a single subset. 
        
        #scoring_function_list = [ConstraintsHandler.e_density_scaled_dev,ConstraintsHandler.chi,ConstraintsHandler.prob_weighted_stat]
        #scoring_function_list = [ConstraintsHandler.prob_weighted_stat,ConstraintsHandler.prob_weighted_stat,ConstraintsHandler.e_density_scaled_dev,ConstraintsHandler.e_density_scaled_dev]
        #scoring_function_list = [ConstraintsHandler.prob_weighted_stat]
        #scoring_function_list = [ConstraintsHandler.prob_weighted_stat,ConstraintsHandler.chi,ConstraintsHandler.e_density_scaled_dev]
        #scoring_function_list = [ConstraintsHandler.chi]
        #scoring_function_list = [ConstraintsHandler.dev_sqr]
        scoring_function_list = [ConstraintsHandler.log_chi]
        #scoring_function_list = [ConstraintsHandler.prob_weighted_stat]
        #scoring_function_list = [ConstraintsHandler.e_density_scaled_dev]
        scoring_function = scoring_function_list[self.loop%len(scoring_function_list)]
        print(f"Using scoring function: {scoring_function.__name__}")

        if self.debug_skip_holton_data_generation:
            #Override
            need_to_prepare_geom_files=False
            
        # keep altloc_subset None if all subsets. This is to save regenerating some files and to make it clear from file names when we are using a subset.
        num_altlocs = len(altloc_subset) if altloc_subset is not None else len(self.model_protein_altlocs)
        do_water_swaps=(self.model_protein_altlocs==self.model_solvent_altlocs) and num_altlocs==2
        if need_to_prepare_geom_files:
            self.prepare_pdb_and_read_altlocs(model_to_swap,model_to_swap,sep_chain_format=False) 
            Solver.LP_Input.prepare_geom_files(model_to_swap,[altloc_subset],
                                                  waters = True,
                                                  water_swaps=do_water_swaps)
            need_to_prepare_geom_files=False

        swapper.clear_candidates()
        if not read_prior_run:
            atoms, connections = Solver.LP_Input(model_to_swap, restrained_refine_pdb_file_path, tensions, self.symmetries, ignore_waters=False,altloc_subset=altloc_subset,resnums=resnums,resnames=resnames).calculate_paths(
                scoring_function=scoring_function,
                clash_punish_thing=False,
                nonbonds=True,   # Note this won't look at nonbonds with water if ignore_waters=True. 
                water_water_nonbond = not DISABLE_WATER_ALTLOC_OPTIM, 
                constraint_weights=self.weight_factors,
            )
            swaps_file_path = Solver.solve(atoms,connections,out_dir=self.output_dir,
                                            out_handle=self.model_handle,
                                            num_solutions=num_solutions,
                                            force_sulfur_bridge_swap_solutions=False, #True
                                            protein_sites=True, 
                                            inert_water_sites=DISABLE_WATER_ALTLOC_OPTIM,
                                            water_sites= not allot_protein_independent_of_waters,
                                            forbid_altloc_changes=dict(name=forbid_altloc_changes),
                                            forbidden_atom_bond_changes=dict(name=forbidden_atom_bond_changes),
                                            forbidden_atom_any_connection_changes=dict(name=forbidden_atom_any_connection_changes),
                                            MAIN_CHAIN_ONLY=MAIN_CHAIN_ONLY,SIDE_CHAIN_ONLY=SIDE_CHAIN_ONLY,NO_CB_CHANGES=NO_CB_CHANGES,
                                            max_bond_changes=self.max_bond_changes,
                                            NO_O_BOND_CHANGES=(self.loop%self.O_bond_change_period!=0)
                                            #water_sites=False,
                                            )
        else:
            swaps_file_path = Solver.swaps_file_path(out_dir=self.output_dir,out_handle=self.model_handle)
        
        # Translate candidate solutions from LinearOptimizer into swaps lists
        # Try proposing each solution until one is accepted or we run out.
        swapper.add_candidates(swaps_file_path) #
    
        candidate_models:list[str]=[]
        candidate_swaps:list[Swapper.SwapGroup]=[]
        candidate_model_dir = f"{self.output_dir}/{self.model_handle}_swapOptions_{self.loop}/"
        os.makedirs(candidate_model_dir,exist_ok=True)
        
        if num_sols_already_saved_this_loop == 0:
            # Remove files from any previous call
            for file in os.listdir(candidate_model_dir):
                path = candidate_model_dir+file
                assert  os.path.abspath(path) != os.path.abspath(model_to_swap)
                os.remove(path) 

        i = num_sols_already_saved_this_loop
        print(f"Applying swaps to {model_to_swap}")
        while swapper.solutions_remaining()>0: 
            ### Swap on unrestrained model ###
            working_model, swapGroup = swapper.run(model_to_swap)
            
            swap_sequence = [swapGroup]
            if not DISABLE_WATER_ALTLOC_OPTIM and allot_protein_independent_of_waters:
                ## Swap waters with protein altlocs fixed
                atoms, connections = Solver.LP_Input(working_model, tensions, self.symmetries,ignore_waters=False).calculate_paths(
                scoring_function=scoring_function,
                clash_punish_thing=False,
                nonbonds=True,
                constraint_weights=self.weight_factors,
                )
                print("Allotting waters")
                #TODO args in dict so don't repeat massive list of args.
                waters_swapped_path = Solver.solve(atoms,connections,out_dir=self.output_dir,
                                                        out_handle=self.model_handle,
                                                        num_solutions=1,
                                                        force_sulfur_bridge_swap_solutions=False,
                                                        protein_sites=True,
                                                        water_sites=True,
                                                        inert_protein_sites=True, # NOTE
                                                        forbid_altloc_changes=dict(name=forbid_altloc_changes),
                                                        forbidden_atom_bond_changes=dict(name=forbidden_atom_bond_changes),
                                                        forbidden_atom_any_connection_changes=dict(name=forbidden_atom_any_connection_changes),
                                                        MAIN_CHAIN_ONLY=MAIN_CHAIN_ONLY,SIDE_CHAIN_ONLY=SIDE_CHAIN_ONLY,NO_CB_CHANGES=NO_CB_CHANGES,
                                                        max_bond_changes=self.max_bond_changes,
                                                        NO_O_BOND_CHANGES=(self.loop%self.O_bond_change_period!=0)
                                                        )
                water_swapper = Swapper()
                water_swapper.add_candidates(waters_swapped_path)
                waterSwapGroup = []
                if len(water_swapper.swap_groups_sorted[0].swaps)!=0:
                    print("Found better water altloc allotments")
                    working_model, waterSwapGroup = water_swapper.run(working_model)
                else:
                    print("Waters unchanged")
                swap_sequence.append(waterSwapGroup)
            
                
            moved_path =  f"{candidate_model_dir }{i+1}.pdb"
            candidate_models.append(moved_path)
            candidate_swaps.append(swap_sequence)
            shutil.move(working_model,moved_path) 
            i+=1
        return candidate_models,candidate_swaps
    
    def determine_best_model(self,model_dir, altloc_subsets_list=None, minimize_R=True,minimize_wE=True,regenerate_R=False):
        assert minimize_R or minimize_wE # both is okay too

        models = os.listdir(model_dir)
        models = [f"{model_dir}{m}"  for m in models]
        best_both_decrease = np.inf
        best = np.inf
        best_tie_breaker_meas=np.inf

        global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
        
        models_for_scoring=[]
        for i, model in enumerate(models):
            altloc_subset = None if altloc_subsets_list is None else altloc_subsets_list[i] 
            if altloc_subset is None:
                modified_model = models[i]
            else:
                struct=PDBParser().get_structure("struct",models[i])
                ordered_atom_lookup = OrderedAtomLookup(struct.get_atoms(),
                                                            protein=True,waters=True,
                                                            altloc_subset=altloc_subset)   
                #temp_path=working_model[:-4]+"_subsetOut.pdb"
                temp_path=Solver.LP_Input.subset_model_path(models[i],altloc_subset)[:-4]+"Out.pdb"
                ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=models[i],out_path=temp_path)
                modified_model=temp_path
            models_for_scoring.append(modified_model)


        def pooled_method(i):
            create_score_file(models_for_scoring[i],self.output_dir,ignore_H=PROPOSE_IGNORES_H,
                              reflections_for_R = self.hkl_path if regenerate_R else None)

        with Pool(self.num_threads) as p:
            p.map(pooled_method,range(len(models)))
                

        best_model_idx=None
        for i, model_for_scoring in enumerate(models_for_scoring):
            print(model)
            combined, wE, Rwork, Rfree = get_score(score_file_name(model_for_scoring,ignore_H=PROPOSE_IGNORES_H))
            #print("Python read | model score wE Rwork Rfree | ",model,combined, wE, Rwork, Rfree)
            if minimize_wE != minimize_R:
                meas, tie_breaker_meas = (Rfree, wE) if minimize_R else (wE,Rfree) 
                if ((meas < best)
                or (meas == best and tie_breaker_meas < best_tie_breaker_meas)
                ):
                    best = meas
                    best_tie_breaker_meas=tie_breaker_meas
                    best_model_idx = i


            else: 
                # Find best Combined Score, but where both decrease relative to current score if possible.

                if Rfree < self.current_score.R_free and wE < self.current_score.wE:
                    assert combined < self.current_score.combined
                    if combined < best_both_decrease:
                        best_both_decrease = combined
                        best_model_idx = i
                if best_both_decrease == np.inf and combined < best:
                    best = combined # this score is no longer used if a result is found where both wE and Rfree decrease
                    best_model_idx = i
        best_model = models[best_model_idx]
        best_untangler_score=Untangler.Score(*get_score(score_file_name(models_for_scoring[i],ignore_H=PROPOSE_IGNORES_H),verbose=False))
        print(f"Best: {best_model} ({best_untangler_score})")
        return best_model
                



    class Strategy(Enum):
        Batch=0
        SwapManyPairs=1
        



    def refinement_loop(self,two_swaps=False,allot_protein_independent_of_waters=False,strategy=None): 
        # working_model stores the path of the current model. It changes value during the loop.
        # TODO if stuck (tried all candidate swap sets for the loop), do random flips or engage Metr.Hastings or track all the new model scores and choose the best.
        

        # if strategy is None:
        #     if len(self.model_protein_altlocs) <= 3:
        #         strategy=Untangler.Strategy.Batch
        #     else:
        #         strategy=Untangler.Strategy.SwapManyPairs
        strategy=Untangler.Strategy.Batch
        #strategy=Untangler.Strategy.SwapManyPairs

        working_model = self.current_model
        if not self.never_do_unrestrained:        
            skip_unrestrained = self.debug_skip_refine or self.debug_skip_unrestrained_refine or (self.loop==self.first_loop and self.debug_skip_first_unrestrained_refine)
            #if (self.loop+9)%10==0 and not skip_unrestrained: # XXX
            if (self.loop+1)%10==0 and not skip_unrestrained: # XXX
                self.working_model = self.nice_long_refine(working_model)

            working_model = self.refine_for_positions(working_model,debug_skip=skip_unrestrained) 
            if skip_unrestrained and not os.path.exists(working_model):
                working_model = self.current_model

                
        num_best_solutions=min(self.loop+self.n_best_swap_start,self.n_best_swap_max) # increase num solutions we search over time...
        
        assert working_model[-4:]==".pdb"
        old_working_model=working_model
        working_model = old_working_model[:-4]+"_fmtd.pdb"
        self.prepare_pdb_and_read_altlocs(old_working_model,working_model)

        def get_tensions(restrained_model,unrestrained_model):
            if not (skip_unrestrained and os.path.exists(geo_file_name(unrestrained_model))): # XXX risky..
                create_score_file(unrestrained_model,log_out_folder_path=self.output_dir)
            if not os.path.exists(geo_file_name(restrained_model)): 
                create_score_file(restrained_model,log_out_folder_path=self.output_dir)
            tensions = GeoXrayTension([restrained_model,unrestrained_model],self.symmetries,water_water_nonbond=not DISABLE_WATER_ALTLOC_OPTIM).site_tensions

            print("Highest tension sites:")
            highest_sites = sorted(tensions, key=tensions.get, reverse=True)[:10]
            for key in highest_sites:
                print(key, tensions[key])

            plt.hist(list(tensions.values()),bins=50)
            plt.savefig("tmp.png")
            plt.close()


            # Negative tensions are set to 0
            for key in tensions:
                tensions[key]=max(0,tensions[key])
            ###
            return tensions

        measure_preswap_postswap = False
        score_file_needs_generation=True
        if measure_preswap_postswap:
            preswap_score = Untangler.Score(*assess_geometry_wE(working_model,self.output_dir))
        if strategy == Untangler.Strategy.Batch:
            
            ## 
            altloc_subsets=[]
            kwargs_list = []


            # Limited options
            if (self.loop+2)%3!=0: # All options
                subset_size=self.altloc_subset_size
                num_combinations=1
                altloc_subsets.extend(self.get_altloc_subsets(subset_size,num_combinations))
                kwargs_list.extend([{} for _ in altloc_subsets])

            else: # Focused swaps V2
                #focused_subset_size =7 
                focused_subset_size =self.altloc_subset_size 
                num_focused_combinations=2
                focused_subsets = self.get_altloc_subsets(focused_subset_size,num_focused_combinations)
                for subset in focused_subsets:
                    #for key in ["NO_CB_CHANGES","MAIN_CHAIN_ONLY"]:
                    for key in ["MAIN_CHAIN_ONLY"]:
                            print(f"Considering {key}")
                            altloc_subsets.append(subset)
                            kwargs_list.append({key:True})


            



            # if self.loop==6:
            #      altloc_subsets=[["B","A","C"]]
            cand_models,cand_swaps = [],[]
            num_sols_stored=0  # FIXME only will load latest file if skipping first swaps. Should load all models in the swapOptions directory.
            last_altloc_subset="NOT SET" # NOTE None means all altlocs!
            tensions=None
            altloc_subsets_list=[] # tracks the altloc_subset used in generating each candidate model 
            for altloc_subset,kwargs in zip(altloc_subsets,kwargs_list):
                skip_swaps = self.debug_skip_first_swaps and self.loop==self.first_loop
                ## TENSIONS ##  #TODO make the one generated for post unrestrained be the same path as generated by candidate_models_From_Swapper so don't need to reuse.
                def get_subset_model(full_model):
                    struct=PDBParser().get_structure("struct",full_model)
                    ordered_atom_lookup = OrderedAtomLookup(struct.get_atoms(),
                                                                protein=True,
                                                                waters=True,
                                                                altloc_subset=altloc_subset)   
                    #temp_path=working_model[:-4]+"_subsetOut.pdb"
                    out_path=Solver.LP_Input.subset_model_path(full_model,altloc_subset)[:-4]+"Out.pdb"
                    ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=full_model,out_path=out_path)
                    return out_path
                need_to_prepare_geom_models=True 
                if not skip_swaps:
                    if altloc_subset!=last_altloc_subset: 
                        #TODO uncomment out once get_subset_model generates subset model to same path as candidate_models_from_swapper AND the model generated is the same (e.g. if excluding water, it should too).
                        unrestrained_subset_model = get_subset_model(self.current_model); #need_to_prepare_geom_models=False
                        restrained_subset_model = get_subset_model(working_model)
                        if TENSIONS:
                            tensions = get_tensions(unrestrained_subset_model,restrained_subset_model)
                            assert tensions is not None, "??"
                    if TENSIONS:
                        assert tensions is not None,(altloc_subset,last_altloc_subset)
                else:
                    need_to_prepare_geom_models=False 
                ###############

                loop_cand_models,loop_cand_swaps = self.candidate_models_from_swapper(self.swapper,num_best_solutions,tensions,working_model,self.current_model, allot_protein_independent_of_waters,
                                        need_to_prepare_geom_files=need_to_prepare_geom_models,read_prior_run=skip_swaps,
                                        altloc_subset=altloc_subset,num_sols_already_saved_this_loop=num_sols_stored,**kwargs)
                altloc_subsets_list.extend((altloc_subset,)*len(loop_cand_models))
                cand_models.extend(loop_cand_models)
                cand_swaps.extend(loop_cand_swaps)
                num_sols_stored+=len(loop_cand_models)
                last_altloc_subset=altloc_subset


            # TODO refinements should be initiated after solution found.

            refined_model_dir = self.regular_batch_refine(cand_models,altloc_subsets_list=altloc_subsets_list,
                                                          debug_skip=self.debug_skip_refine)
            #working_model = self.determine_best_model(refined_model_dir)
            working_model = self.determine_best_model(refined_model_dir,altloc_subsets_list=altloc_subsets_list, 
                                                      minimize_wE=CONSIDER_WE_WHEN_CHOOSING_BEST_BATCH,
                                                      regenerate_R=False)
            
            #### TODO Sucks make better ####
            best_model_that_was_refined = os.path.basename(working_model).split("_")[-1]
            candidate_model_dir = f"{self.output_dir}/{self.model_handle}_swapOptions_{self.loop}/"
            best_model_that_was_refined = candidate_model_dir+best_model_that_was_refined
            swaps = cand_swaps[cand_models.index(best_model_that_was_refined)]
            ################################
            
            print("Refining all conformers for best model")
            if altloc_subsets is None or all([v is None for v in altloc_subsets]):
                shutil.copy(working_model,self.get_out_path(f"loopEnd{self.loop}"))
            else:
                working_model=self.regular_refine(working_model,num_loops_override=1, debug_skip=self.debug_skip_refine)
            score_file_needs_generation=False
            if measure_preswap_postswap:
                postswap_score = Untangler.Score(*assess_geometry_wE(best_model_that_was_refined,self.output_dir))
                print("Score preswap:",preswap_score) 
                print("Score postswap:",postswap_score) 
        elif strategy == Untangler.Strategy.SwapManyPairs:
            tensions = get_tensions(self.current_model,working_model) if TENSIONS else None

            swaps_focused = None
            # FOCUS SWAPS (TODO better name)
            if (not self.debug_skip_first_focus_swaps) or self.loop!=self.first_loop:
                #working_model,swaps_cys = self.swap_cysteines(self.swapper,working_model)
                working_model,swaps_cys = self.swap_cysteines(self.swapper,tensions,working_model)
                #working_model,swaps_tyr = self.swap_specific_residues(["TYR"],self.swapper,tensions,working_model,file_tag="TyrSwaps")
                swaps_focused = swaps_cys#+swaps_tyr

            ##

            independent_sulfur_approach=False
            forbidden_atom_any_connection_changes=[]
            if independent_sulfur_approach:
                forbidden_atom_any_connection_changes.append("SG")
            working_model,swaps = self.many_swapped(self.swapper,tensions,working_model,self.current_model,allot_protein_independent_of_waters,
                                                    forbidden_atom_any_connection_changes=forbidden_atom_any_connection_changes,
                                                    altloc_subset_size=self.altloc_subset_size)
            if swaps_focused is not None:
                swaps = swaps_focused+swaps
            if measure_preswap_postswap:
                postswap_score = Untangler.Score(*assess_geometry_wE(working_model,self.output_dir))
                print("Score preswap:",preswap_score) 
                print("Score postswap:",postswap_score) 
            working_model=self.regular_refine(working_model,debug_skip=self.debug_skip_refine)
        else:
            raise Exception(f"Invalid strategy {strategy}")
        
        new_model_was_accepted = self.propose_model(working_model,score_file_needs_generation=score_file_needs_generation)
        
        '''
        if not new_model_was_accepted and two_swaps:
            #### Swap again on restrained model. Take Best solution only.  ###
            unswapper=Swapper()
            cand_models,cand_swaps = self.candidate_models_from_swapper(unswapper,1,working_model,allot_protein_independent_of_waters)

            raise Exception("Unimplemented")
            if "unswapper.swap_groups_sorted[0].swaps has a non water swap": # TODO or if unswaps == swaps
                working_model = self.determine_best_model(refined_model_dir)
                swaps += cand_swaps[cand_models.index(working_model)]
                new_model_was_accepted = self.propose_model(working_model)

            else:
                print("Continuing, no unswaps found")
        '''

        if new_model_was_accepted:
            self.swaps_history.append(swaps)
            self.write_swaps_history()
        
    def write_swaps_history(self):
        with open(f"{self.output_dir}refine_logs/swapHistory_{self.model_handle}.txt","w") as f:
            f.writelines('\n'.join([str(swap) for swap in self.swaps_history]))

    def propose_model(self,working_model,score_file_needs_generation=True):
        if score_file_needs_generation or PROPOSE_IGNORES_H:
            create_score_file(working_model,log_out_folder_path=self.output_dir,ignore_H=PROPOSE_IGNORES_H)
        new_score = Untangler.Score(*get_score(score_file_name(working_model,ignore_H=PROPOSE_IGNORES_H)))
        
        #deltaE = new_wE-self.current_score.wE # geometry only
        if CONSIDER_WE_WHEN_CHOOSING_BEST_BATCH:
            deltaE = new_score.combined-self.current_score.combined  # include R factor
            max_wE_increase = self.max_wE_frac_increase*self.current_score.wE
        else:
            deltaE = new_score.R_free-self.current_score.R_free  # include R factor
            max_wE_increase = self.max_wE_frac_increase*self.current_score.R_free
        if self.current_score.wE==np.inf:
            max_wE_increase = np.inf
        p_accept = self.P_accept_swap(
            deltaE,
            max_wE_increase=max_wE_increase
        )
        outcome,set_new_model="xXx Rejected xXx",False
        if random.random() < p_accept or self.debug_always_accept_proposed_model:
            self.current_model = working_model
            #accepted_path=self.get_out_path(f"Accepted{self.loop}")
            self.current_score = new_score
            if self.best_score.combined > self.current_score.combined:
                print(f"previous_best: {self.best_score}")
                self.best_score = self.current_score
            outcome,set_new_model="oOo Accepted oOo",True
            shutil.copy(working_model,self.get_out_path(f"Accepted{self.loop}"))
        print(f"{outcome} proposed model change with score of {new_score} (P_accept: {p_accept:.2f}) ")
        UntangleFunctions.copy_model_and_geo(self.current_model,self.next_current_model_name())
        self.current_model = self.next_current_model_name()
        return set_new_model



    def current_model_name_at_loop(self,i):
        return self.get_out_path(f"current{i}")
    def next_current_model_name(self):
        return self.current_model_name_at_loop(self.loop+1)
        return Untangler.output_dir+f"{self.model_handle}_current_{self.loop}.pdb"

    def initial_refine(self,model_path,**kwargs)->str:
        # Try to get atoms as close to their true positions as possible
        for wc, wu, n_cycles in zip([1,0.5,0.2,0.1],[1,0,0,0],[2,4,5,5]):
        #for wc, wu, n_cycles in zip([1,0.5],[1,0],[8,4]):
        #for wc, wu, n_cycles in zip([1],[1],[self.num_end_loop_refine_cycles]):
            if self.refinement==self.PHENIX:
                refine_params = self.get_refine_params_phenix(
                    "initial",
                    model_path=model_path,
                    num_macro_cycles=n_cycles,
                    wc=wc,
                    wu=wu,
                    hold_water_positions=self.holding_water(),
                    refine_occupancies=False,
                    ordered_solvent=PHENIX_ORDERED_SOLVENT or self.debug_phenix_ordered_solvent_on_initial,
                    shake=0,
                    refine_hydrogens=True, # XXX
                )
            elif self.refinement == self.REFMAC:
                refine_params=self.get_refine_params_refmac(
                    "initial",
                    model_path=model_path,
                    unrestrained=False,
                    # max_trials=100,
                    max_trials=150,
                    min_trials=0,
                    dampA=0.1,
                    dampB=0.1,
                    dampC=0.25,
                    wc=0.5,
                    refine_water_occupancies=self.refmac_refine_water_occupancies_initial
                )
            model_path = self.refine(
                refine_params,
                **kwargs
            )
        return model_path

    def refine_for_positions(self,model_path,**kwargs)->str:

        # TODO CRITICAL measure how geo changes. Use this to weight geo in altloc assignment (if it changes little, weight it less).

        next_model=model_path
        # No idea why need to do this loop. But otherwise it jumps at 1_xyzrec
        if self.refinement==self.PHENIX:
            for n in range(self.num_refine_for_positions_macro_cycles_phenix):
                refine_params = self.get_refine_params_phenix(
                    #f"unrestrained-mc{self.loop}-{n}",
                    f"unrestrained{self.loop}",
                    model_path=next_model,
                    num_macro_cycles=1,
                    #wc=0,
                    wc=self.refine_for_positions_geo_weight,
                    max_sigma_movement_waters=0.1,#0.001, 
                    #wu=0,
                    #wc=0.25,
                    hold_water_positions=True,
                    refine_hydrogens=True, # XXX
                )
                next_model = self.refine(
                    refine_params,
                    **kwargs
                )
        elif self.refinement == self.REFMAC:
            refine_params=self.get_refine_params_refmac(
                f"unrestrained{self.loop}",
                model_path=next_model,
                #unrestrained=True,
                unrestrained=False,
                # dampA=0.004,
                # dampB=0.01,
                # dampA=0.02,
                # dampB=0.04,
                wc=1,
                # dampA=0.5,
                # dampB=0.5,
                # dampC=0.5,
                dampA=0.05,
                dampB=0.25,
                dampC=0.5,
                min_trials=1,
                max_trials=3,
            )
            next_model = self.refine(
                refine_params,
                **kwargs
            )

        return next_model
    
        
    #TODO regular_batch_refine and regular_refine should get refine params from same source 
    def regular_batch_refine(self,model_paths:list[str],altloc_subsets_list=None, **kwargs):
        param_set = []
        for i, model in enumerate(model_paths):
            altloc_subset = None if (altloc_subsets_list is None) else altloc_subsets_list[i]
            if self.refinement==self.PHENIX:
                refine_params = self.get_refine_params_phenix(
                    f"loopEnd{self.loop}-{i+1}",
                    model_path=model,
                    num_macro_cycles=self.num_end_loop_refine_cycles,
                    wc=self.wc_anneal_start if self.wc_anneal_loops==0 else min(1,self.wc_anneal_start+(self.loop/self.wc_anneal_loops)*(1-self.wc_anneal_start)),
                    hold_water_positions=self.holding_water(),
                    #refine_occupancies=False,
                    #ordered_solvent=False,
                    ordered_solvent=PHENIX_ORDERED_SOLVENT,
                    refine_occupancies=False,
                    #max_sigma_movement_waters=0.05,
                    max_sigma_movement_waters=0.1,
                    altloc_subset=altloc_subset,
                    #max_sigma_movement_waters=0.07,
                    )
            elif self.refinement == self.REFMAC:
                refine_params=self.get_refine_params_refmac(
                    f"loopEnd{self.loop}-{i+1}",
                    model_path=model,
                    unrestrained=False,
                    refine_water_occupancies=True,
                    min_trials=1,
                    max_trials=10,
                    # dampA=0.1,
                    # dampB=0.25,
                    # dampA=1,
                    # dampB=1,
                    # dampA=0.6,
                    # dampB=0.6,
                    # dampC=0.8,
                    dampA=0.05,
                    dampB=0.25,
                    dampC=0.5,
                    wc=0.5,
                    altloc_subset=altloc_subset
                )
            param_set.append(refine_params)
        return self.batch_refine(f"loopEnd{self.loop}",param_set,**kwargs)


    def nice_long_refine(self,model_path,**kwargs)->str:
        print("Performing nice long restrained refinement ")
        if self.refinement==self.PHENIX:
            refine_params_H_refine = self.get_refine_params_phenix(
                f"loopStart{self.loop}",
                model_path=model_path,
                num_macro_cycles=2 * self.num_end_loop_refine_cycles,
                refine_hydrogens=True, # NOTE  the hydrogens get rotated around to the wrong place, so need to do this every now and then. Subsequent riding H refine idealizes the hydrogen ADP and brings Rfree back up. (Possible to get to solution without for challenge without doing this, but you may have to get lucky and also have to do enough refines for the hydrogens to slowly rotate around to the right places)
                hold_water_positions=True
            )
            working_model = self.refine(
                refine_params_H_refine,
                **kwargs
            )
            refine_params = self.get_refine_params_phenix(
                f"loopStart{self.loop}",
                model_path=working_model,
                num_macro_cycles=2 * self.num_end_loop_refine_cycles,
                refine_hydrogens=False, 
            )
            return self.refine(
                refine_params,
                **kwargs
            )
        elif self.refinement == self.REFMAC:
            print("Not implemented for refmac")
            return model_path


    def regular_refine(self,model_path,altloc_subset=None,num_loops_override=None,**kwargs)->str:
        print("Performing post-reallotment restrained refinement ")
        if self.refinement==self.PHENIX:
            refine_params = self.get_refine_params_phenix(
                f"loopEnd{self.loop}",
                model_path=model_path,
                num_macro_cycles=self.num_end_loop_refine_cycles if num_loops_override is None else num_loops_override,
                wc= self.wc_anneal_start if self.wc_anneal_loops==0 else min(self.default_wc,self.wc_anneal_start+(self.loop/self.wc_anneal_loops)*(self.default_wc-self.wc_anneal_start)),
                hold_water_positions=self.holding_water(),
                #refine_occupancies=False,
                #ordered_solvent=False,
                ordered_solvent=PHENIX_ORDERED_SOLVENT,
                refine_occupancies=False,
                #max_sigma_movement_waters=0.05,
                max_sigma_movement_waters=0.1,
                altloc_subset=altloc_subset,
                #max_sigma_movement_waters=0.07,

            )
        elif self.refinement == self.REFMAC:
            old_model_path = model_path
            model_path = old_model_path[:-4]+"_sepChn.pdb"
            self.prepare_pdb_and_read_altlocs(old_model_path,model_path,sep_chain_format=True)

            refine_params=self.get_refine_params_refmac(
                f"loopEnd{self.loop}",
                model_path=model_path,
                unrestrained=False,
                refine_water_occupancies=True,
                min_trials=1,
                max_trials=30,
                # dampA=0.1,
                # dampB=0.25,
                # dampA=1,
                # dampB=1,
                # dampA=0.6,
                # dampB=0.6,
                # dampC=0.8,
                dampA=0.05,
                dampB=0.25,
                dampC=0.5,
                wc=0.5
            )
        return self.refine(
            refine_params,
            **kwargs
        )

    def holding_water(self)->bool:
        return self.loop<self.num_loops_water_held
    def get_out_path(self,out_tag):
        return f"{Untangler.output_dir}{self.model_handle}_{out_tag}.pdb"
    def refine(self,refine_params:(tuple[SimpleNamespace,list[str]]),debug_skip=False,show_python_params=False)->str:
        # assert model_path[-4:]==".pdb", model_path
        P, args = refine_params
        out_path = self.get_out_path(P.out_tag)
        assert os.path.exists(refine_params[0].model_path)
        if not debug_skip:
            max_attempts=3
            attempt=0
            backup_path = out_path+'#'
            if os.path.abspath(refine_params[0].model_path)!=os.path.abspath(out_path): # So as not to disrupt multiple refines using same file output name
                if os.path.exists(out_path):
                    shutil.move(out_path,backup_path)
            while True:
                if show_python_params:
                    print (f"Params: {P}")
                print (f"|+ Running: {' '.join(args)}")
                subprocess.run(args)#,stdout=log)

                if os.path.exists(out_path): #TODO replace with direct way to check for success
                    break
                elif attempt < max_attempts:
                    attempt+=1
                    print(f"Warning: refinement failed for unknown reason! Retrying...")
                    sleep(2)
                else:
                    raise Exception(f"refinement failed {max_attempts} times!")
        
            # NOTE if somethin goes wrong here, the out_path is a file copied from Refinement/tmp_refinement/, so it can be recovered.
            self.prepare_pdb_and_read_altlocs(out_path,out_path+"tmp",sep_chain_format=False) 
            shutil.move(out_path+"tmp",out_path)

        
        return out_path
    
    def group_waters_to_altlocs(self,model_path,out_path):
        print("Grouping waters")
        k=len(self.model_protein_altlocs)
        
        structure:Structure = PDBParser(model_path).get_structure("struct",model_path)
        ordered_ordered_waters:list[Atom]=[]
        water_coords=[]
        water_residues={}
        starting_resnum=np.inf
        for atom in structure.get_atoms():
            if res_is_water(atom.get_parent()):
                resnum=OrderedAtomLookup.atom_res_seq_num(atom)
                starting_resnum=int(min(starting_resnum,resnum))
                water_residues[resnum]=atom.get_parent()
                for ordered_atom in atom:
                    ordered_atom:Atom
                    ordered_ordered_waters.append(atom)
                    water_coords.append(ordered_atom.get_coord())

        
        neighbours = NearestNeighbors(n_neighbors=k).fit(np.array(water_coords))
        distances, indices = neighbours.kneighbors(water_coords)
        indices:list[list[int]] = indices.tolist()
        print("Warning: water altloc grouping code is not implemented correctly")
        for g, group in enumerate(indices):
            assert len(group)==len(self.model_protein_altlocs)
            for i, altloc in enumerate(self.model_protein_altlocs):
                ordered_ordered_waters[group[i]].set_altloc(altloc)
                ordered_ordered_waters[group[i]].set_occupancy(1/k)
                ordered_ordered_waters[group[i]].set_parent(water_residues[int(g+starting_resnum)])

        io = PDBIO()
        io.set_structure(structure)
        io.save(out_path)             

        self.model_solvent_altlocs=self.model_protein_altlocs


    def get_refine_params_refmac(self,out_tag,model_path,unrestrained=False,refine_water_occupancies=False, turn_off_bulk_solvent=TURN_OFF_BULK_SOLVENT,
                                 altloc_subset=None,
                                 min_trials=0,max_trials=None,
                                 wc=0.5,dampA=0.01,dampB=0.05,dampC=0.1):
        assert altloc_subset is None, "altloc subset not implemented!"
        
        reflections_path = self.hkl_path
        param_dict = locals()
        del param_dict["self"]
        P = SimpleNamespace(**param_dict)

        args=["bash", 
            f"{self.working_dir}/{self.refine_refmac_shell_file}",f"{P.model_path}",f"{P.reflections_path}",
            "-o",f"{self.model_handle}_{P.out_tag}",
        ]
        if P.min_trials > 0:
            args.extend(["-m",f"{P.min_trials}"])
        if P.max_trials is not None:
            args.extend(["-n",f"{P.max_trials}"])
        args.extend([
            "-A",f"{P.dampA}",
            "-B",f"{P.dampB}",
            "-C",f"{P.dampC}",
            "-c",f"{P.wc}",
        ])
        for bool_param, flag in ([P.unrestrained,"u"],[P.refine_water_occupancies,"W"],[P.turn_off_bulk_solvent,"t"]):
            if bool_param:
                args.append(f"-{flag}")
        return P, args


            
    def get_refine_params_phenix(self, out_tag=None, model_path=None,num_macro_cycles=None, # mandatory
                          wc=None,wu=1., shake=0., optimize_R=False,
                          hold_water_positions=False,hold_protein_positions=False,
                          refine_occupancies=False,turn_off_bulk_solvent=TURN_OFF_BULK_SOLVENT,ordered_solvent=False,
                          no_restrain_movement=False,max_sigma_movement_waters=0.1,refine_hydrogens=False, # restraining movement refers to the reference_coordinate_restraints option
                          altloc_subset=None):
        if altloc_subset is not None:
            altloc_subset = ''.join(altloc_subset)  
        if wc is None:
            wc = self.default_wc
        ### Override next_model with formatted one.
        #next_model = model_path[:-4]+"_fmtd.pdb"

        # next_model=model_path[:-4]+"_sep_chains.pdb"
        # # Create formatted
        # self.prepare_pdb_and_read_altlocs(model_path,next_model,sep_chain_format=True)

        #model_path = next_model 
        ###
        if PHENIX_FREEZE_WATER:
            hold_water_positions=True
            
        disable_CDL = PHENIX_DISABLE_CDL

        assert (not hold_water_positions) or (not hold_protein_positions) 
        
        reflections_path = self.hkl_path
        param_dict = locals()
        del param_dict["self"]
        #del param_dict["next_model"]
        #####################

        P = SimpleNamespace(**param_dict)
  
        args=["bash", 
            f"{self.working_dir}/{self.refine_shell_file}",f"{P.model_path}",f"{P.reflections_path}",
            "-c",f"{P.wc}",
            "-u",f"{P.wu}",
            "-n",f"{P.num_macro_cycles}",
            "-o",f"{self.model_handle}_{P.out_tag}",
            "-s",f"{P.shake}",
            "-q",f"{P.max_sigma_movement_waters}",
        ]
        if altloc_subset is not None:
            args+= ["-a",altloc_subset]
        for bool_param, flag in ([P.hold_water_positions,"-h"],[P.refine_hydrogens,"-H"],[P.optimize_R,"-r"],
                                 [P.hold_protein_positions,"-p"],[P.refine_occupancies,"-O"],[P.turn_off_bulk_solvent,"-t"],
                                 [P.ordered_solvent,"-S"],[P.no_restrain_movement,"-R"],[P.disable_CDL,"-C"]):
            if bool_param:
                args.append(flag)
        return P, args

    def batch_refine(self,batch_tag,refine_arg_sets:list[dict[str]],debug_skip=False)->str:
        out_directory = f"{self.output_dir}/{self.model_handle}_{batch_tag}/"
        os.makedirs(out_directory,exist_ok=True) 
        # Remove files from any previous call
        if not debug_skip:
            for file in os.listdir(out_directory):
                os.remove(out_directory + file) 

        # TODO folder in output/refine_logs/
        
            global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
            def pooled_method(i):
                max_attempts=3
                attempt=0
                while True:
                    sleep(0.2*i) # Desperate attempt to reduce phenix seg faults.
                    print(f"Refining {i+1}/{len(refine_arg_sets)}")
                    out_path = self.refine(refine_arg_sets[i])
                    out_directory = f"{self.output_dir}/{self.model_handle}_{batch_tag}/"
                    if os.path.exists(out_path):
                        shutil.move(out_path,f"{out_directory}{batch_tag}_{i+1}.pdb") 
                        break
                    elif attempt < max_attempts:
                        attempt+=1
                        print(f"Warning: refinement {i} failed for unknown reason! Retrying...")
                    else:
                        print(f"Warning: refinement {i} failed {max_attempts} times! Skipping...")
                        break
            with Pool(self.num_threads) as p:
                p.map(pooled_method,range(len(refine_arg_sets)))
        return out_directory




    def P_accept_swap(self,deltaE,max_wE_increase,metropolis_hastings=True): 
        if deltaE > max_wE_increase:
            return 0
        if not metropolis_hastings: # allow small (max_wE_increase) increases
            return 1
        k=1 # TODO find good default
        return np.exp(-deltaE/(k*self.acceptance_temperature))


    def refine_for_positions_costly(self,model_path,**kwargs)->str:
        #TODO reduce num cycles if Rfree doesn't decrease.
        #  

        refine_param_sets  = []


        #model_path_A=model_path_B=model_path_C=model_path_D=model_path_E=model_path

        refine_param_sets.append(dict(
            #num_macro_cycles=3,
            num_macro_cycles=2,
            wc=0,
            hold_water_positions=True,
            **kwargs
        ))

        assert False, "TODO: Sequential refine in batch"
        refine_param_sets.append(
            [dict(
            num_macro_cycles=2,
            wc=0.1,
            shake=0.03,
            hold_water_positions=True,
            **kwargs
            ),

            dict(
            num_macro_cycles=1,
            wc=0,
            wu=1,
            hold_water_positions=True,
            **kwargs)
            ]
        )

        # model_path_C = self.refine(
        #     model_path,
        #     "unrestrained_C",
        #     num_macro_cycles=2,
        #     optimize_R=True,
        #     hold_water_positions=True,
        #     **kwargs
        # )
        refine_param_sets.append(dict(
            model_path,
            "unrestrained_D",
            num_macro_cycles=2,
            optimize_R=True,
            hold_water_positions=self.holding_water(),
            **kwargs
        ))

        # model_path_E = self.refine(
        #     model_path,
        #     "unrestrained_E",
        #     num_macro_cycles=3,
        #     wc=0,
        #     shake=0.03,
        #     hold_water_positions=True,
        #     **kwargs
        # )
        refine_param_sets.append(dict(
            model_path,
            "unrestrained_F",
            #num_macro_cycles=5,
            num_macro_cycles=3,
            wc=0.1,
            wu=0,
            hold_water_positions=self.holding_water(),
            **kwargs
        ))
        out_dir = self.batch_refine(model_path,"unrestrained_options",refine_param_sets)



        lowest_Rfree_model=None
        lowest_Rfree = np.inf
        R_dict = {}
        reflections_path=pdb_data_dir()+"refme.mtz" # NOTE
        # with open(f"{self.output_dir}refine_logs/{self.model_handle}_position_refine_Rs.txt","w") as f:
        #     #for m_path in (model_path_A,model_path_B,model_path_C,model_path_D,model_path_E,model_path_F):
        #     for m_path in (os.listdir(out_dir)):
        #         Rwork,Rfree = get_R(m_path,reflections_path)
        #         R_dict[m_path] = f"Rwork: {Rwork}, Rfree: {Rfree}"
        #         if Rfree < lowest_Rfree:
        #             lowest_Rfree_model=m_path
        #             lowest_Rfree=Rfree
        #             print(f"Current positions model (Rwork | Rfree) - ({Rwork} | {Rfree})")
        #         f.write(f"{m_path.split('/')[-1]} | {R_dict[m_path]}\n")
        best_model = self.determine_best_model(out_dir,minimize_wE=False)

        return lowest_Rfree_model




def main():
    if len(sys.argv)!=3:
        print("Usage: python3.9 untangle.py data/myInitialModel.pdb data/myReflections.mtz")
        return

    starting_model = sys.argv[1]
    xray_data = sys.argv[2]
    Untangler(
        # max_num_best_swaps_considered=5,
        default_wc=1,
        num_end_loop_refine_cycles=6,
        max_num_best_swaps_considered=15,
        starting_num_best_swaps_considered=15,
        altloc_subset_size=2,
        #refine_for_positions_geo_weight=0.03,
        refine_for_positions_geo_weight=0,
        num_refine_for_positions_macro_cycles_phenix=1,
        # weight_factors = {
        #     ConstraintsHandler.BondConstraint: 150,
        #     ConstraintsHandler.AngleConstraint: 0.25,
        #     #ConstraintsHandler.NonbondConstraint: 1,  # TODO experiment with this.
        #     ConstraintsHandler.NonbondConstraint: 0,  # TODO experiment with this.
        #     ConstraintsHandler.ClashConstraint: 10,
        # },
        #max_bond_changes=2,
        max_bond_changes=9999,
        # weight_factors = {
        #     ConstraintsHandler.BondConstraint: 1,
        #     ConstraintsHandler.AngleConstraint: 1, #80
        #     ConstraintsHandler.NonbondConstraint: 0,  # TODO experiment with this.
        #     ConstraintsHandler.ClashConstraint: 1,
        # },
        # weight_factors = {
        #     ConstraintsHandler.BondConstraint: 1,
        #     ConstraintsHandler.AngleConstraint: 80,
        #     ConstraintsHandler.NonbondConstraint: 1,  # TODO experiment with this.
        #     ConstraintsHandler.ClashConstraint: 1,
        # },
        # weight_factors = {
        #     ConstraintsHandler.BondConstraint: 1,
        #     ConstraintsHandler.AngleConstraint: 1,
        #     ConstraintsHandler.NonbondConstraint: 1,  # TODO experiment with this.
        #     ConstraintsHandler.ClashConstraint: 0,
        # },
        weight_factors = {
            ConstraintsHandler.BondConstraint: 0.1,
            ConstraintsHandler.AngleConstraint: 80,#1,
            ConstraintsHandler.NonbondConstraint: 0.1,  # TODO experiment with this.
            # ConstraintsHandler.ClashConstraint: 0.1,
            # ConstraintsHandler.TwoAtomPenalty: 1,
            ConstraintsHandler.ClashConstraint: 0,  #0 1 100
            ConstraintsHandler.TwoAtomPenalty: 0,
        },
        # weight_factors = {
        #     ConstraintsHandler.BondConstraint: 1,
        #     ConstraintsHandler.AngleConstraint: 100,
        #     ConstraintsHandler.NonbondConstraint: 0,  # TODO experiment with this.
        #     ConstraintsHandler.ClashConstraint: 0,
        # },
        ).run(
        starting_model,
        xray_data,
        desired_score=18.59,
        max_num_runs=100,
    )
if __name__=="__main__":
    main()

# %%
