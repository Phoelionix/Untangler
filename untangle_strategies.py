#%%
from typing import Any
from LinearOptimizer import Solver
from LinearOptimizer.Input import OrderedAtomLookup
from LinearOptimizer.Swapper import Swapper
import UntangleFunctions
from UntangleFunctions import assess_geometry_wE, get_R,pdb_data_dir,create_score_file,get_score,score_file_name,res_is_water, parse_symmetries_from_pdb
import subprocess
import os, sys
import numpy as np
import shutil
import random
from types import SimpleNamespace 
from multiprocessing import Pool, Process
from time import sleep
from enum import Enum
import itertools
from Bio.PDB import PDBParser,Structure,PDBIO
from Bio.PDB.Atom import Atom,DisorderedAtom
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KNeighborsClassifier
from untangle import Untangler


class UntangleStrategy():
    
    def __init__(self,untangler:Untangler):
        self.untangler = untangler

    def found_better()->str:
        assert False, "Unimplemented"

    # Keep searching for next-best solutions until we see an improvement
    def FindSolutionsAndRefineOnDemand(max_loops=50):
        max_queue_size=12
        process_queue:list[Process] = []
        found_sol = None
        for l in range(max_loops):
            
            args
            process_queue.append(
                Process(target=refine_method, args=args)
            )
            process_queue[-1].start()
            

            
            if len(process_queue)>max_queue_size:
                process_queue[0].join()
                del process_queue[0]

            found_sol = found_better(out_model_dir)
            if found_sol is not None:
                break
        
        for p in process_queue:
            p.join()
        del process_queue
        
        if found_sol is None:
            # Check the remaining processes
            found_sol = found_better(out_model_dir)
        
        return found_sol

        ###


# class ParamPresets():


# %%
