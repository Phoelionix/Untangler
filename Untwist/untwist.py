#%%
import os, sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
import shutil
from Untwist.detect_twists import create_untwist_file
from Untwist.apply_untwists import apply_untwists 
from Bio.PDB import PDBParser,Structure
from Bio.PDB.Atom import DisorderedAtom, Atom
from LinearOptimizer.Tag import DisorderedTag
import numpy as np

#TODO
# Should focus on finding solutions that are as close to 90 degrees as possible. 


#(I) they share the same midpoint
# (II) they are farther apart
# (III) each conformer position is related via rotation of a model conformer around the vector formed by a and b.
# 
# 
# 
# 

def apply_all_untwists_for_two_conformations(working_model):
    untwist_file = create_untwist_file(working_model)
    untwist_move_models_dir, all_untwists_model, changes_only_model = apply_untwists(working_model,untwist_file)

    return all_untwists_model, changes_only_model

def get_untwist_atom_options(working_model)->list[DisorderedAtom]:
    untwist_file = create_untwist_file(working_model)
    untwist_move_models_dir, all_untwists_model, changes_only_model = apply_untwists(working_model,untwist_file)
    #TODO perform unrestrained refinement with only the atoms changed allowed to move. See which "stick".
    disordered_atoms = list(PDBParser().get_structure("struct",changes_only_model).get_atoms())
    assert all(type(v)==DisorderedAtom for v in disordered_atoms)
    return disordered_atoms

def get_untwist_atom_options_that_survived_unrestrained(pos_refined_model, pre_untwist_model, changes_only_model,max_gap_close_frac,exclude_H):
    untwist_atoms = {DisorderedTag.from_atom(a):a for a in PDBParser().get_structure("struct",changes_only_model).get_atoms()}
    original_atoms = {DisorderedTag.from_atom(a):a for a in PDBParser().get_structure("struct",pre_untwist_model).get_atoms()
                      if DisorderedTag.from_atom(a) in untwist_atoms}
    post_refine_atoms = {DisorderedTag.from_atom(a):a for a in PDBParser().get_structure("struct",pos_refined_model).get_atoms()
                         if DisorderedTag.from_atom(a) in untwist_atoms}
    

    def separation(a:Atom,b:Atom):
        return np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
    
    alternate_atoms:list[DisorderedAtom]=[]
    disallowed:list[DisorderedAtom]=[]
    for site_tag , untwist_atm in untwist_atoms.items():
        passed=True
        if exclude_H and untwist_atm.element=="H":
            continue
        og_atm = original_atoms[site_tag]
        post_ref_atm=post_refine_atoms[site_tag]

        # Find closest atoms before unrestrained refine.
        # Very lazy check TODO
        pre_refined_min_sep=np.inf
        post_refined_min_sep=np.inf
        for ordered_og_atm in og_atm:
            for ordered_untwist_atm in untwist_atm:
                pre_refined_min_sep = min(pre_refined_min_sep,separation(ordered_og_atm,ordered_untwist_atm))
            for ordered_post_ref_atm in post_ref_atm:
                post_refined_min_sep = min(post_refined_min_sep,separation(ordered_og_atm,ordered_post_ref_atm))

        frac_gap_closed = abs(post_refined_min_sep-pre_refined_min_sep)/pre_refined_min_sep
        if frac_gap_closed>max_gap_close_frac:
            passed = False
        if separation(*post_ref_atm)<separation(*og_atm):
            passed = False
        if passed:
            alternate_atoms.append(post_ref_atm)
        else:
            disallowed.append(post_ref_atm)
                
    return alternate_atoms,disallowed
        



    



    
# %%
