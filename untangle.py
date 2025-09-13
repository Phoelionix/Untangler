#%%
from typing import Any
from LinearOptimizer import Solver
from LinearOptimizer.Input import OrderedAtomLookup
from LinearOptimizer.Swapper import Swapper
from UntangleFunctions import assess_geometry_wE, get_R,pdb_data_dir,create_score_file,get_score,score_file_name,res_is_water
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


class Untangler():
    working_dir = os.path.abspath(os.getcwd())
    output_dir = f"{working_dir}/output/"
    refine_shell_file=f"Refinement/Refine.sh"
    ####
    # Skip stages. Requires files to already have been generated (up to the point you are debugging).
    debug_skip_refine = False
    debug_skip_initial_refine=True
    debug_skip_first_unrestrained_refine=True
    debug_skip_first_swaps=False
    debug_skip_unrestrained_refine=False
    debug_skip_holton_data_generation=False
    debug_skip_initial_holton_data_generation=debug_skip_initial_refine
    ####
    num_threads=10
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
                 wc_anneal_start=1,wc_anneal_loops=0, starting_num_best_swaps_considered=20, # 20,
                 max_num_best_swaps_considered=100,num_unrestrained_macro_cycles=3,
                 num_loops_water_held=0):
        self.set_hyper_params(acceptance_temperature,max_wE_frac_increase,num_end_loop_refine_cycles,
                              wc_anneal_start,wc_anneal_loops, starting_num_best_swaps_considered,
                              max_num_best_swaps_considered,num_unrestrained_macro_cycles, 
                              num_loops_water_held)
        self.previously_swapped = []
        self.model_handle=None
        self.current_model=None
        self.initial_score=self.current_score=self.best_score=None
        self.best_score=self.current_score
        self.swapper:Swapper=None
        self.swaps_history:list[Swapper.SwapGroup] = [] 
        self.model_protein_altlocs=None
        self.model_solvent_altlocs=None
    def set_hyper_params(self,acceptance_temperature=1,max_wE_frac_increase=0, num_end_loop_refine_cycles=2,  # 8,
                 wc_anneal_start=1,wc_anneal_loops=0, starting_num_best_swaps_considered=5, # 20,
                 max_num_best_swaps_considered=100,num_unrestrained_macro_cycles=3,
                 num_loops_water_held=0):
        # NOTE Currently max wE increase is percentage based, 
        # but TODO optimal method needs to be investigated.
        self.n_best_swap_start=starting_num_best_swaps_considered
        self.n_best_swap_max = max_num_best_swaps_considered
        self.num_unrestrained_macro_cycles=num_unrestrained_macro_cycles
        self.max_wE_frac_increase=max_wE_frac_increase
        self.acceptance_temperature=acceptance_temperature
        self.num_end_loop_refine_cycles=num_end_loop_refine_cycles
        self.num_loops_water_held=num_loops_water_held
        self.wc_anneal_start = wc_anneal_start
        self.wc_anneal_loops=wc_anneal_loops

    def prepare_pdb_and_read_altlocs(self,pdb_path,out_path,sep_chain_format=False):
        # Gets into format we expect. !!!!!!Assumes single chain!!!!!
        protein_altlocs = []
        solvent_altlocs = []
        with open(pdb_path) as I:
            start_lines = []
            end_lines = []
            atom_dict:dict[str,dict[str,dict[str,str]]] = {}  
            last_chain=None
            solvent_res_names=["HOH"]
            solvent_chain_id = "z"
            warned_collapse=False
            def replace_chain(line,chain_id):
                chain_id = str(chain_id)
                assert len(chain_id)==1
                return line[:21]+chain_id+line[22:]
            def replace_serial_num(line,serial_num):
                serial_num = str(serial_num)
                serial_num = ' '*(5-len(serial_num))+serial_num
                return line[:6]+serial_num+line[11:]

            for line in I:
                if line.startswith("TER"):
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
                resnum = line[22:26]
                if resname in solvent_res_names:
                    end_lines+=replace_chain(line,solvent_chain_id)
                    if altloc not in solvent_altlocs:
                        solvent_altlocs.append(altloc) 
                    continue
                assert len(end_lines)==0
                    
                if not sep_chain_format and not warned_collapse and chain != last_chain and last_chain is not None:
                    print("Warning: Multiple chains detected. Collapsing chains into single chain")
                    warned_collapse=True
                if resnum not in atom_dict:
                    atom_dict[resnum] = {}
                if altloc not in atom_dict[resnum]:
                    atom_dict[resnum][altloc] = {}
                    
                    if altloc not in protein_altlocs:
                        protein_altlocs.append(altloc) 
                         
                atom_dict[resnum][altloc][name]=line  
                last_chain = chain
                continue
                    
            n=0
            if not sep_chain_format: # format for untangler stuff
                protein_chain_id = "A"
                for res_atom_dict in atom_dict.values():
                    for altloc_atom_dict in res_atom_dict.values():
                        for line in altloc_atom_dict.values():
                            n+=1
                            modified_line = replace_chain(line,protein_chain_id)
                            modified_line = replace_serial_num(modified_line,n)
                            start_lines.append(modified_line)
            else: # Note that lines for each chain need to be contiguous in the file
                chain_dict={}
                for res_atom_dict in atom_dict.values():
                    for altloc, altloc_atom_dict in res_atom_dict.items():
                        protein_chain_id = altloc
                        for line in altloc_atom_dict.values():
                            n+=1
                            modified_line = replace_chain(line,protein_chain_id)
                            modified_line = replace_serial_num(modified_line,n)
                            if altloc not in chain_dict:
                                chain_dict[altloc]=[]
                            chain_dict[altloc].append(modified_line)
                for _, lines in chain_dict.items():
                    for modified_line in lines:
                        start_lines.append(modified_line)

        with open(out_path,'w') as O:
            O.writelines(start_lines+end_lines)
        self.model_protein_altlocs=protein_altlocs
        self.model_solvent_altlocs=solvent_altlocs

    def run(self,pdb_file_path,hkl_file_path,desired_score=18.6,max_num_runs=100):
        # pdb_file_path: path to starting model
        # hkl_file_path: path to reflection data.
        # TODO Currently assume in data folder.
        assert hkl_file_path[-4:]==".mtz", f"hkl path doesn't end in '.mtz': {hkl_file_path}"
        self.hkl_handle = os.path.basename(hkl_file_path)[:-4] 
        #assert os.path.dirname(hkl_file_path)[-5:-1]=="data", hkl_file_path    
        assert pdb_file_path[-4:]==".pdb", f"file path doesn't end in '.pdb': {pdb_file_path}"
        self.model_handle = os.path.basename(pdb_file_path)[:-4]
        self.current_model= Untangler.output_dir+f"{self.model_handle}_current.pdb"
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
        if self.model_solvent_altlocs != self.model_protein_altlocs:
            self.group_waters_to_altlocs(self.current_model,self.current_model)
            assert self.model_solvent_altlocs==self.model_protein_altlocs
        self.swapper = Swapper()


        self.loop=0
        initial_model=self.initial_refine(self.current_model,debug_skip=(self.debug_skip_refine or self.debug_skip_initial_refine)) 
        if not self.debug_skip_holton_data_generation and not self.debug_skip_initial_holton_data_generation:
            self.initial_score = Untangler.Score(*assess_geometry_wE(initial_model, log_out_folder_path=self.output_dir))
        else:
            self.initial_score = Untangler.Score(*get_score(score_file_name(initial_model)))
        self.current_score = self.best_score = Untangler.Score.inf_bad_score()# i.e. force accept first solution
        shutil.copy(initial_model,self.current_model)
        
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


    def many_swapped(self,swapper,model_to_swap:str,allot_protein_independent_of_waters:bool,altloc_subset_size=2,num_combinations=15,repeats=2):
        #TODO try strategy of making one altloc as good as possible, while other can be terrible.
        
        working_model = f"{self.output_dir}/{self.model_handle}_manySwaps.pdb"
        
        all_swaps=[]
        if self.debug_skip_first_swaps and self.loop == 0:
            return working_model, ["Unknown"]
        
        shutil.copy(model_to_swap,working_model)
        if len(self.model_protein_altlocs)<=altloc_subset_size:
            repeats = 1 # No point repeating if considering all at once.
        for r in range(repeats+1):

            altloc_subset_combinations:itertools.combinations[tuple[str]] = list(itertools.combinations(self.model_protein_altlocs, altloc_subset_size))
            print("num possible combinations:",len(altloc_subset_combinations))
            if len(altloc_subset_combinations) > num_combinations:
                altloc_subset_combinations = random.sample(altloc_subset_combinations,num_combinations)
            for n, altloc_subset in enumerate(altloc_subset_combinations): 
                header=f"=========Altloc Allotment {n+1}/{len(altloc_subset_combinations)}, Cycle {r+1}/{repeats+1}=========="
                print(f"\n{header}\n{'-'*len(header)}")
                print(f"Optimizing connections across altlocs {', '.join(altloc_subset)}")
                cand_models, cand_swaps = self.candidate_models_from_swapper(swapper,1,working_model,allot_protein_independent_of_waters,
                                                                   altloc_subset=altloc_subset)
                assert len(cand_models)==len(cand_swaps)==1
                shutil.move(cand_models[0],working_model)
                all_swaps.extend(cand_swaps[0])

                measure_wE_after=False
                if measure_wE_after:
                    struct=PDBParser().get_structure("struct",working_model)
                    ordered_atom_lookup = OrderedAtomLookup(struct.get_atoms(),
                                                                protein=True,waters=not allot_protein_independent_of_waters,
                                                                altloc_subset=altloc_subset)   
                    temp_path=working_model[:-4]+"_subsetOut.pdb"
                    ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=working_model,out_path=temp_path)
                    assess_geometry_wE(temp_path,self.output_dir)

        print("=======End Altloc Allotment=============\n")
        return working_model, all_swaps
    def candidate_models_from_swapper(self,swapper:Swapper,num_solutions,model_to_swap:str,allot_protein_independent_of_waters:bool,altloc_subset=None,skip_geom_file_generation=False): #TODO refactor as method of Swapper class
        # TODO should be running solver for altloc set partitioned into subsets, not a single subset. 
        
        if self.debug_skip_holton_data_generation:
            #Override
            skip_geom_file_generation=True
        swapper.clear_candidates()
        atoms, connections = Solver.MTSP_Solver(model_to_swap,ignore_waters=allot_protein_independent_of_waters,altloc_subset=altloc_subset,skip_subset_file_generation=skip_geom_file_generation).calculate_paths(
            clash_punish_thing=False,
            nonbonds=True,   # Note this won't look at nonbonds with water if ignore_waters=True. 
            water_water_nonbond = False,  # This is so we don't get a huge number of nearly identical solutions from swapping waters around.
            skip_geom_file_generation=skip_geom_file_generation,
        )
        swaps_file_path = Solver.solve(atoms,connections,out_dir=self.output_dir,
                                        out_handle=self.model_handle,
                                        num_solutions=num_solutions,
                                        force_sulfur_bridge_swap_solutions=False, #True
                                        protein_sites=True, 
                                        water_sites= not allot_protein_independent_of_waters,
                                        #water_sites=False,
                                        )
        
        # Translate candidate solutions from LinearOptimizer into swaps lists
        # Try proposing each solution until one is accepted or we run out.
        swapper.add_candidates(swaps_file_path) #
    
        candidate_models:list[str]=[]
        candidate_swaps:list[Swapper.SwapGroup]=[]
        candidate_model_dir = f"{self.output_dir}/{self.model_handle}_swapOptions_{self.loop}/"
        os.makedirs(candidate_model_dir,exist_ok=True)
        # Remove files from any previous call
        for file in os.listdir(candidate_model_dir):
            path = candidate_model_dir+file
            assert  os.path.abspath(path) != os.path.abspath(model_to_swap)
            os.remove(path) 

        i = 0
        while swapper.solutions_remaining()>0: 
            ### Swap on unrestrained model ###
            working_model, swapGroup = swapper.run(model_to_swap)
            
            swap_sequence = [swapGroup]
            if allot_protein_independent_of_waters:
                ## Swap waters with protein altlocs fixed
                atoms, connections = Solver.MTSP_Solver(working_model,ignore_waters=False).calculate_paths(
                clash_punish_thing=False,
                nonbonds=True 
                )
                print("Allotting waters")
                waters_swapped_path = Solver.solve(atoms,connections,out_dir=self.output_dir,
                                                        out_handle=self.model_handle,
                                                        num_solutions=1,
                                                        force_sulfur_bridge_swap_solutions=False,
                                                        protein_sites=True,
                                                        water_sites=True,
                                                        inert_protein_sites=True, # NOTE
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
    
    def determine_best_model(self,model_dir, minimize_R=True,minimize_wE=True):
        assert minimize_R or minimize_wE # both is okay too

        models = os.listdir(model_dir)
        models = [f"{model_dir}{m}"  for m in models]
        best_both_decrease = np.inf
        best = np.inf
        best_model=None

        global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
        def pooled_method(i):
            create_score_file(self.output_dir,models[i])

        with Pool(self.num_threads) as p:
            p.map(pooled_method,range(len(models)))
                

        for model in models:
            print(model)
            combined, wE,Rwork, Rfree = get_score(score_file_name(model))
            print("Python read | model score wE Rwork Rfree | ",model,combined, wE, Rwork, Rfree)
            if minimize_wE != minimize_R:
                meas = Rfree if minimize_R else wE 
                if meas < best:
                    best = meas
                    best_model = model

            else: 
                # Find best Combined Score, but where both decrease relative to current score if possible.

                if Rfree < self.current_score.R_free and wE < self.current_score.wE:
                    assert combined < self.current_score.combined
                    if combined < best_both_decrease:
                        best_both_decrease = combined
                        best_model = model
                if best_both_decrease == np.inf and combined < best:
                    best = combined # this score is no longer used if a result is found where both wE and Rfree decrease
                    best_model = model
        assert best_model is not None
        print("Best:",best_model)
        return best_model
                



    class Strategy(Enum):
        Batch=0
        SwapManyPairs=1
        



    def refinement_loop(self,two_swaps=False,allot_protein_independent_of_waters=False,strategy=None): 
        # working_model stores the path of the current model. It changes value during the loop.
        # TODO if stuck (tried all candidate swap sets for the loop), do random flips or engage Metr.Hastings or track all the new model scores and choose the best.
        
        if strategy is None:
            if len(self.model_protein_altlocs) == 2:
                strategy=Untangler.Strategy.Batch
            else:
                strategy=Untangler.Strategy.SwapManyPairs
                
        working_model = self.refine_for_positions(self.current_model,debug_skip=self.debug_skip_refine or self.debug_skip_unrestrained_refine or (self.loop==0 and self.debug_skip_first_unrestrained_refine)) 
        
            
        num_best_solutions=min(self.loop+self.n_best_swap_start,self.n_best_swap_max) # increase num solutions we search over time...
        
        self.prepare_pdb_and_read_altlocs(working_model,working_model)
        preswap_score = Untangler.Score(*assess_geometry_wE(working_model,self.output_dir))
        if strategy == Untangler.Strategy.Batch:
            if self.debug_skip_first_swaps:
                assert False, "Unimplemented"
            cand_models,cand_swaps = self.candidate_models_from_swapper(self.swapper,num_best_solutions,working_model,allot_protein_independent_of_waters)
            refined_model_dir = self.regular_batch_refine(cand_models,debug_skip=self.debug_skip_refine)
            working_model = self.determine_best_model(refined_model_dir)
            #### TODO Sucks make better ####
            best_model_that_was_refined = os.path.basename(working_model).split("_")[-1]
            candidate_model_dir = f"{self.output_dir}/{self.model_handle}_swapOptions_{self.loop}/"
            best_model_that_was_refined = candidate_model_dir+best_model_that_was_refined
            ################################
            postswap_score = Untangler.Score(*assess_geometry_wE(best_model_that_was_refined,self.output_dir))
            print("Score preswap:",preswap_score) 
            print("Score postswap:",postswap_score) 
            swaps = cand_swaps[cand_models.index(best_model_that_was_refined)]
        elif strategy == Untangler.Strategy.SwapManyPairs:
            working_model,swaps = self.many_swapped(self.swapper,working_model,allot_protein_independent_of_waters)
            postswap_score = Untangler.Score(*assess_geometry_wE(working_model,self.output_dir))
            print("Score preswap:",preswap_score) 
            print("Score postswap:",postswap_score) 
            working_model=self.regular_refine(working_model,debug_skip=self.debug_skip_refine)
        else:
            raise Exception(f"Invalid strategy {strategy}")
        
        new_model_was_accepted = self.propose_model(working_model)
        
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

    def propose_model(self,working_model):
        
        new_score = Untangler.Score(*assess_geometry_wE(working_model,self.output_dir))
        #deltaE = new_wE-self.current_score.wE # geometry only
        deltaE = new_score.combined-self.current_score.combined  # include R factor
        max_wE_increase = self.max_wE_frac_increase*self.current_score.wE
        if self.current_score.wE==np.inf:
            max_wE_increase = np.inf
        p_accept = self.P_accept_swap(
            deltaE,
            max_wE_increase=max_wE_increase
        )
        outcome,set_new_model="xXx Rejected xXx",False
        if random.random() < p_accept:
            self.current_model = self.next_current_model_name()
            shutil.copy(working_model,self.current_model)
            self.current_score = new_score
            if self.best_score.combined > self.current_score.combined:
                print(f"previous_best: {self.best_score}")
                self.best_score = self.current_score
            outcome,set_new_model="oOo Accepted oOo",True
        print(f"{outcome} proposed model change with score of {new_score} (P_accept: {p_accept:.2f}) ")
        return set_new_model

    def next_current_model_name(self):
        return Untangler.output_dir+f"{self.model_handle}_current_{self.loop}.pdb"

    def initial_refine(self,model_path,**kwargs)->str:
        # Try to get atoms as close to their true positions as possible
        #for wc, wu, n_cycles in zip([1,0.5,0.2,0.1],[1,0,0,0],[2,4,5,5]):
        #for wc, wu, n_cycles in zip([1,0.5],[1,0],[8,4]):
        for wc, wu, n_cycles in zip([1],[1],[self.num_end_loop_refine_cycles]):
            refine_params = self.get_refine_params(
                "initial",
                model_path=model_path,
                num_macro_cycles=n_cycles,
                wc=wc,
                wu=wu,
                hold_water_positions=self.holding_water(),
                #ordered_solvent=True
            )
            model_path = self.refine(
                refine_params,
                **kwargs
            )
        return model_path

    def refine_for_positions(self,model_path,**kwargs)->str:
        # No idea why need to do this. But otherwise it jumps at 1_xyzrec
        next_model=model_path
        for n in range(self.num_unrestrained_macro_cycles):
            refine_params = self.get_refine_params(
                f"unrestrained-mc{n}",
                model_path=next_model,
                num_macro_cycles=1,
                wc=0,
                hold_water_positions=True,
            )
            next_model = self.refine(
                refine_params,
                **kwargs
            )
        return next_model
    
        
    #TODO regular_batch_refine and regular_refine should get refine params from same source 
    def regular_batch_refine(self,model_paths,**kwargs):
        param_set = []
        for i, model in enumerate(model_paths):
            param_set.append(self.get_refine_params(
                f"loopEnd{self.loop}-{i+1}",
                model_path=model,
                num_macro_cycles=self.num_end_loop_refine_cycles,
                wc=self.wc_anneal_start if self.wc_anneal_loops==0 else min(1,self.wc_anneal_start+(self.loop/self.wc_anneal_loops)*(1-self.wc_anneal_start)),
                hold_water_positions=self.holding_water(),
                #ordered_solvent=True,
                #refine_occupancies=False
                ))
        return self.batch_refine(f"loopEnd{self.loop}",param_set,**kwargs)


    def regular_refine(self,model_path,**kwargs)->str:

        refine_params = self.get_refine_params(
            f"loopEnd{self.loop}",
            model_path=model_path,
            num_macro_cycles=self.num_end_loop_refine_cycles,
            wc= self.wc_anneal_start if self.wc_anneal_loops==0 else min(1,self.wc_anneal_start+(self.loop/self.wc_anneal_loops)*(1-self.wc_anneal_start)),
            hold_water_positions=self.holding_water(),
            #ordered_solvent=True,
            #refine_occupancies=False
        )
        return self.refine(
            refine_params,
            **kwargs
        )
        # return self.refine(
        #     model_path,
        #     "loopEnd",
        #     num_macro_cycles=self.num_end_loop_refine_cycles,
        #     wc=min(1,self.wc_anneal_start+(self.loop/self.wc_anneal_loops)*(1-self.wc_anneal_start)),
        #     #hold_water_positions= self.holding_water(),  
        #     hold_protein_positions = True,
        #     **kwargs
        # )

    def holding_water(self)->bool:
        return self.loop<self.num_loops_water_held
    def refine(self,refine_params:(tuple[SimpleNamespace,list[str]]),debug_skip=False,show_python_params=False)->str:
        # assert model_path[-4:]==".pdb", model_path
        P, args = refine_params
        out_path = f"{self.output_dir}/{self.model_handle}_{P.out_tag}.pdb"
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
                print (f"Running: {' '.join(args)}")
                subprocess.run(args)#,stdout=log)

                if os.path.exists(out_path): #TODO replace with direct way to check for success
                    break
                elif attempt < max_attempts:
                    attempt+=1
                    print(f"Warning: refinement failed for unknown reason! Retrying...")
                else:
                    raise Exception(f"refinement failed {max_attempts} times!")
                

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



    #def get_refine_params(self,param_dict: SimpleNamespace | dict[str,Any]):
    def get_refine_params(self, out_tag=None, model_path=None,num_macro_cycles=None, # mandatory
                          wc=1,wu=1, shake=0, optimize_R=False,
                          hold_water_positions=False,hold_protein_positions=False,
                          refine_occupancies=True,turn_off_bulk_solvent=True,ordered_solvent=False):
        ### Override next_model with formatted one.
        #next_model = model_path[:-4]+"_fmtd.pdb"
        #self.prepare_pdb_and_read_altlocs(model_path,next_model,sep_chain_format=True)
        #model_path = next_model 
        ###

        assert (not hold_water_positions) or (not hold_protein_positions) 
        param_dict = locals()
        del param_dict["self"]
        #del param_dict["next_model"]
        #####################


        #TODO remove.
        P_defaults = dict(out_tag=None,model_path=None,num_macro_cycles=None, # mandatory
                          wc=1,wu=1, shake=0, optimize_R=False,
                          hold_water_positions=False,hold_protein_positions=False,
                          refine_occupancies=True,turn_off_bulk_solvent=True,ordered_solvent=False)
        
        if type(param_dict) is SimpleNamespace:
            param_dict:dict[str,Any] = vars(param_dict)
        for key, default_value in P_defaults.items():
            if key not in param_dict:
                param_dict[key] = default_value 
        assert P_defaults.keys() == param_dict.keys() 
        P = SimpleNamespace(**param_dict)
  
        args=["bash", 
            f"{self.working_dir}/{self.refine_shell_file}",f"{P.model_path}",f"{P.out_tag}",
            "-c",f"{P.wc}",
            "-u",f"{P.wu}",
            "-n",f"{P.num_macro_cycles}",
            "-o",f"{self.model_handle}_{P.out_tag}",
            "-s",f"{P.shake}",
            "-d", self.hkl_handle]
        for bool_param, flag in ([P.hold_water_positions,"-h"],[P.optimize_R,"-r"],[P.hold_protein_positions,"-p"],[P.refine_occupancies,"-O"],[P.turn_off_bulk_solvent,"-t"],[P.ordered_solvent,"-S"]):
            if bool_param:
                args.append(flag)
        return P, args

    def batch_refine(self,batch_tag,refine_arg_sets:list[dict[str]],debug_skip=False)->str:
        out_directory = f"{self.output_dir}/{self.model_handle}_{batch_tag}/"
        os.makedirs(out_directory,exist_ok=True) 
        # Remove files from any previous call
        for file in os.listdir(out_directory):
            os.remove(out_directory + file) 

        # TODO folder in output/refine_logs/
        
        global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
        def pooled_method(i):
            max_attempts=3
            attempt=0
            while True:
                sleep(2*i) # Desperate attempt to reduce phenix seg faults.
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
        if not debug_skip:
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
        max_num_best_swaps_considered=30,
        starting_num_best_swaps_considered=10,
        num_unrestrained_macro_cycles=1
        ).run(
        starting_model,
        xray_data,
        desired_score=18.41,
        max_num_runs=100
    )
if __name__=="__main__":
    main()

# %%
