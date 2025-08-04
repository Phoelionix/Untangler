import sys 
#%%
from Bio.PDB import PDBParser#, parse_pdb_header
from Bio.PDB.Structure import Structure#, parse_pdb_header
from Bio.PDB.Atom import Atom,DisorderedAtom#, parse_pdb_header
from Bio.PDB.MMCIFParser import MMCIFParser#, parse_pdb_header
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
import sys
from typing import List
from copy import deepcopy
import os
import subprocess
from residue_grouper import ResidueGrouper
from InputFromModel import InputFromModel
import random
import numpy as np

def refine_worst_residues(model_handle,from_one_end=False):
    model_handle

    # pdbPath = "/home/speno/Untangle/StructureGeneration/data/" \
    # f"{handle}.pdb"
    holton_worst_geo_file = "/path/to/HoltonOutputs/" \
    f"{model_handle}_worstatoms.txt"

    num_worst=6
    if from_one_end:
        should_ignore_last_residue_refined #if little change when focusing specific refine of that residue then assume it's correct.
        num_worst=10
    bad_atoms = InputFromModel.HoltonWorstGeo(holton_worst_geo_file,
                num_worst=num_worst,min_badness=7).bad_atoms
    
    bad_atoms_bonds = InputFromModel.HoltonWorstGeo(holton_worst_geo_file,
                num_worst=num_worst,min_badness=2,badness_types=["bond"]).bad_atoms
    bad_atoms+=bad_atoms_bonds
    bad_atoms = set(bad_atoms)
    first_only=False 
    bad_residues = []
    for name,res,badness,kind in bad_atoms:
        if res in bad_residues:
            continue
        bad_residues.append(res)
        if first_only:
            break

    if len(bad_residues)==0:
        return
    #e.g. "resid 43 or resid 44 or resid 35"
    refinement_selection = f"resid {bad_residues[0]}"
    for n in bad_residues[1:]:
        refinement_selection+=f" or resid {n}"
    print(refinement_selection)
    refine_shell_file ="/path/to/RefineSpecific.sh"
    # wc=random.uniform(0.01,2)
    # wu=random.uniform(0.01,2)


    args=["bash", f"{refine_shell_file}",f"{model_handle}",refinement_selection,"-o",f"{model_handle}_RefWrst"]
    print (f"$Running {args}")
    subprocess.run(args)#,stdout=log)
    print("finished")
    print("--==--")
if __name__=="__main__":
    refine_worst_residues(sys.argv[1])

# %%
