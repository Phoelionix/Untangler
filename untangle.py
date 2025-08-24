#%%
from LinearOptimizer import Solver
from LinearOptimizer.Swapper import Swapper
from UntangleFunctions import assess_geometry_wE, get_R,pdb_data_dir
import subprocess
import os, sys
import numpy as np
import shutil
import random



class Untangler():
    working_dir = os.path.abspath(os.getcwd())
    output_dir = f"{working_dir}/output/"
    refine_shell_file=f"Refinement/Refine.sh"
    ####
    debug_skip_refine = False # Skip refinement stages. Requires files to already have been generated (up to the point you are debugging).
    debug_skip_initial_refine=False 
    class Score():
        def __init__(self,combined,wE,R_work):
            self.wE=wE
            self.combined=combined  
            self.R_work=R_work
        def __repr__(self):
            return f"{self.combined} | wE: {self.wE} Rwork: {self.R_work}"
        @staticmethod
        def inf_bad_score()->'Untangler.Score':
            return Untangler.Score(np.inf,np.inf,np.inf)

    def __init__(self,acceptance_temperature=1,max_wE_frac_increase=0, num_end_loop_refine_cycles=5, 
                 wc_anneal_start=0.6,wc_anneal_loops=2, starting_num_best_swaps_considered=20,
                 max_num_best_swaps_considered=100,num_loops_water_held=1):
        self.set_hyper_params(acceptance_temperature,max_wE_frac_increase,num_end_loop_refine_cycles,
                              wc_anneal_start,wc_anneal_loops, starting_num_best_swaps_considered,
                              max_num_best_swaps_considered,num_loops_water_held)
        self.previously_swapped = []
        self.model_handle=None
        self.current_model=None
        self.initial_score=self.current_score=self.best_score=None
        self.best_score=self.current_score
        self.swapper:Swapper=None
    def set_hyper_params(self,acceptance_temperature=1,max_wE_frac_increase=0, num_end_loop_refine_cycles=5, 
                 wc_anneal_start=0.6,wc_anneal_loops=2, starting_num_best_swaps_considered=20,
                 max_num_best_swaps_considered=100,num_loops_water_held=1):
        # NOTE Currently max wE increase is percentage based, 
        # but TODO optimal method needs to be investigated.
        self.n_best_swap_start=starting_num_best_swaps_considered
        self.n_best_swap_max = max_num_best_swaps_considered
        self.max_wE_frac_increase=max_wE_frac_increase
        self.acceptance_temperature=acceptance_temperature
        self.num_end_loop_refine_cycles=num_end_loop_refine_cycles
        self.num_loops_water_held=num_loops_water_held
        self.wc_anneal_start = wc_anneal_start
        self.wc_anneal_loops=wc_anneal_loops

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
        shutil.copy(pdb_file_path,self.current_model)
        self.swapper = Swapper()

        initial_model=self.initial_refine(self.current_model,debug_skip=(self.debug_skip_refine or self.debug_skip_initial_refine)) 
        self.initial_score = Untangler.Score(*assess_geometry_wE(self.output_dir,initial_model))
        self.current_score = self.best_score = Untangler.Score.inf_bad_score()# i.e. force accept first solution
        shutil.copy(initial_model,self.current_model)
        self.loop=0
        
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
        # working_model = self.refinement_loop()
        # self.propose_model(working_model)
        self.refinement_loop()

    def refinement_loop(self,two_swaps=True):
        # working_model stores the path of the current model. It changes value during the loop.
        # TODO if stuck (tried all candidate swap sets for the loop), do random flips or engage Metr.Hastings or track all the new model scores and choose the best.
        working_model = self.refine_for_positions(self.current_model,debug_skip=self.debug_skip_refine) 
        
        self.swapper.clear_candidates()
        atoms, connections = Solver.MTSP_Solver(working_model, protein_sites=True,water_sites=False).calculate_paths(
            clash_punish_thing=False,
            nonbonds=False   # Use later when swapping waters. i.e. we look for the best protein conformation and allot waters after.
        )
        num_best_solutions=min(self.loop+self.n_best_swap_start,self.n_best_swap_max) # increase num solutions we search over time...
        swaps_file_path = Solver.solve(atoms,connections,out_dir=self.output_dir,
                                        out_handle=self.model_handle,
                                        num_solutions=num_best_solutions,
                                        force_sulfur_bridge_swap_solutions=True)
        
        # Translate candidate solutions from LinearOptimizer into swaps lists
        # Try proposing each solution until one is accepted or we run out.
        self.swapper.add_candidates(swaps_file_path) # Note sorted then added to end. So it will go through one loaded set first, then next if no improvements found there.
    

        pre_swap_model = working_model
        while self.swapper.solutions_remaining()>0: 
            ### Swap on unrestrained model ###
            working_model, _ = self.swapper.run(pre_swap_model)

            # Swap waters
            atoms, connections = Solver.MTSP_Solver(working_model, protein_sites=False,water_sites=True).calculate_paths(
            clash_punish_thing=False,
            nonbonds=True 
            )
            waters_swapped_path = Solver.solve(atoms,connections,out_dir=self.output_dir,
                                                    out_handle=self.model_handle,
                                                    num_solutions=1,
                                                    force_sulfur_bridge_swap_solutions=False)
            water_swapper = Swapper()
            water_swapper.add_candidates(waters_swapped_path)
            if len(water_swapper.swap_groups_sorted[0].swaps)!=0:
                print("Found better water altloc allotments")
                working_model, _ = water_swapper.run(working_model)
            
            working_model=self.regular_refine(working_model,debug_skip=self.debug_skip_refine)
            new_model_was_accepted = self.propose_model(working_model)

            if not new_model_was_accepted and two_swaps:
                #### Swap again on restrained model. Take Best solution only.  ###
                atoms, connections = Solver.MTSP_Solver(working_model,protein_sites=True,water_sites=False).calculate_paths(
                    clash_punish_thing=False,
                    nonbonds=True
                )
        
                unswap_file_path = Solver.solve(atoms,connections,out_dir=self.output_dir,out_handle=self.model_handle,num_solutions=1,force_sulfur_bridge_swap_solutions=False)
                unswapper = Swapper()
                unswapper.add_candidates(unswap_file_path)
                if len(unswapper.swap_groups_sorted[0].swaps)!=0:
                    working_model, _ = unswapper.run(working_model)

                    working_model=self.regular_refine(working_model,debug_skip=self.debug_skip_refine)
                    new_model_was_accepted = self.propose_model(working_model)
                else:
                    print("Continuing, no unswaps found")

            if new_model_was_accepted:
                break
        
        #return working_model 

    def propose_model(self,working_model):
        
        new_combined, new_wE,new_R = assess_geometry_wE(self.output_dir,working_model)
        new_score = Untangler.Score(new_combined,new_wE,new_R)
        #deltaE = new_wE-self.current_score.wE # geometry only
        deltaE = new_combined-self.current_score.combined  # include R factor
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
                self.best_score = self.current_score
            outcome,set_new_model="oOo Accepted oOo",True
        print(f"{outcome} proposed model change with score of {new_score} (P_accept: {p_accept:.2f}) ")
        return set_new_model

    def next_current_model_name(self):
        return Untangler.output_dir+f"{self.model_handle}_current_{self.loop}.pdb"

    def initial_refine(self,model_path,**kwargs)->str:
        # Try to get atoms as close to their true positions as possible
        for wc, wu, n_cycles in zip([1,0.5,0.2,0.1],[1,0,0,0],[2,4,5,5]):
            model_path = self.refine(
                model_path,
                "initial",
                num_macro_cycles=n_cycles,
                wc=wc,
                wu=wu,
                hold_water_positions=True,
                **kwargs
            )
        return model_path

    def refine_for_positions(self,model_path,**kwargs)->str:
        return self.refine(
            model_path,
            "unrestrained",
            num_macro_cycles=1,
            wc=0,
            hold_water_positions=True,
            **kwargs
        )
    
    def regular_refine(self,model_path,**kwargs)->str:
        return self.refine(
            model_path,
            "loopEnd",
            num_macro_cycles=self.num_end_loop_refine_cycles,
            wc=min(1,self.wc_anneal_start+(self.loop/self.wc_anneal_loops)*(1-self.wc_anneal_start)),
            hold_water_positions= self.holding_water(),  # even when False, bond distances are still held (refine_water_bond_length_hold_template.eff)
            **kwargs
        )
    def holding_water(self)->bool:
        return self.loop<self.num_loops_water_held
    def refine(self,model_path,out_tag,num_macro_cycles,wc=1,wu=1,optimize_R=False,hold_water_positions=False,debug_skip=False,shake=0)->str:
        assert model_path[-4:]==".pdb", model_path
        out_model_handle = os.path.basename(model_path)[:-4]
        if not debug_skip:
            args=["bash", 
              f"{self.working_dir}/{self.refine_shell_file}",f"{model_path}",f"{out_tag}",
              "-c",f"{wc}",
              "-u",f"{wu}",
              "-n",f"{num_macro_cycles}",
              "-o",f"{self.model_handle}_{out_tag}",
              "-s",f"{shake}",
              "-d", self.hkl_handle]
            for bool_param, flag in ([hold_water_positions,"-h"],[optimize_R,"-r"]):
                if bool_param:
                    args.append(flag)
            print (f"Running {args}")
            subprocess.run(args)#,stdout=log)
        out_path = f"{self.output_dir}/{self.model_handle}_{out_tag}.pdb"

        return out_path


    def P_accept_swap(self,deltaE,max_wE_increase,metropolis_hastings=True): 
        if deltaE > max_wE_increase:
            return 0
        if not metropolis_hastings: # allow small (max_wE_increase) increases
            return 1
        k=1 # TODO find good default
        return np.exp(-deltaE/(k*self.acceptance_temperature))


    # def refine_for_positions_costly(self,model_path,**kwargs)->str:
    #     #TODO reduce num cycles if Rfree doesn't decrease.
    #     #  
    #     reflections_path=pdb_data_dir()+"refme.mtz" # NOTE

    #     #model_path_A=model_path_B=model_path_C=model_path_D=model_path_E=model_path

    #     #TODO parallelise
    #     model_path_A = self.refine(
    #         model_path,
    #         "unrestrained_A",
    #         #num_macro_cycles=3,
    #         num_macro_cycles=2,
    #         wc=0,
    #         hold_water_positions=True,
    #         **kwargs
    #     )
    #     model_path_B = self.refine(
    #         model_path,
    #         "unrestrained_B",
    #         num_macro_cycles=2,
    #         wc=0.1,
    #         shake=0.03,
    #         hold_water_positions=True,
    #         **kwargs
    #     )
    #     model_path_B = self.refine(
    #         model_path_B,
    #         "unrestrained_B",
    #         num_macro_cycles=1,
    #         wc=0,
    #         wu=1,
    #         hold_water_positions=True,
    #         **kwargs
    #     )

    #     # model_path_C = self.refine(
    #     #     model_path,
    #     #     "unrestrained_C",
    #     #     num_macro_cycles=2,
    #     #     optimize_R=True,
    #     #     hold_water_positions=True,
    #     #     **kwargs
    #     # )
    #     model_path_D = self.refine(
    #         model_path,
    #         "unrestrained_D",
    #         num_macro_cycles=2,
    #         optimize_R=True,
    #         hold_water_positions=self.holding_water(),
    #         **kwargs
    #     )

    #     # model_path_E = self.refine(
    #     #     model_path,
    #     #     "unrestrained_E",
    #     #     num_macro_cycles=3,
    #     #     wc=0,
    #     #     shake=0.03,
    #     #     hold_water_positions=True,
    #     #     **kwargs
    #     # )
    #     model_path_F = self.refine(
    #         model_path,
    #         "unrestrained_F",
    #         #num_macro_cycles=5,
    #         num_macro_cycles=3,
    #         wc=0.1,
    #         wu=0,
    #         hold_water_positions=self.holding_water(),
    #         **kwargs
    #     )
    #     lowest_Rfree_model=None
    #     lowest_Rfree = np.inf
    #     R_dict = {}
    #     with open(f"{self.output_dir}refine_logs/{self.model_handle}_position_refine_Rs.txt","w") as f:
    #         #for m_path in (model_path_A,model_path_B,model_path_C,model_path_D,model_path_E,model_path_F):
    #         for m_path in (model_path_A,model_path_B,model_path_D,model_path_F):
    #             Rwork,Rfree = get_R(m_path,reflections_path)
    #             R_dict[m_path] = f"Rwork: {Rwork}, Rfree: {Rfree}"
    #             if Rfree < lowest_Rfree:
    #                 lowest_Rfree_model=m_path
    #                 lowest_Rfree=Rfree
    #                 print(f"Current positions model (Rwork | Rfree) - ({Rwork} | {Rfree})")
    #             f.write(f"{m_path.split('/')[-1]} | {R_dict[m_path]}\n")
    #     return lowest_Rfree_model




def main():
    if len(sys.argv)!=3:
        print("Usage: python3.9 untangle.py data/myInitialModel.pdb data/myReflections.mtz")
        return

    starting_model = sys.argv[1]
    xray_data = sys.argv[2]
    Untangler().run(
        starting_model,
        xray_data,
        desired_score=18.4,
        max_num_runs=100
    )
if __name__=="__main__":
    main()

# %%
