#%%
import sys,os,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from ConformationTree.split_conformation_geoms import create_all_child_restraints
from LinearOptimizer.Tag import *
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler
from LinearOptimizer.Input import LP_Input
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from Bio.PDB import PDBParser,Structure,PDBIO
from UntangleFunctions import parse_symmetries_from_pdb
from ConformationTree.split_pdb import split_specific


def run(model_path,child_parent_altlocs_dict):
    out_dir = os.path.join(UntangleFunctions.UNTANGLER_WORKING_DIRECTORY,"ConformationTree","output")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    parent_altlocs,_=UntangleFunctions.get_altlocs_from_pdb(model_path)
    parent_altlocs=sorted(list(parent_altlocs))

    excluded_resnames=["CYS",]
    ordered_atom_lookup = OrderedAtomLookup(model_path, waters=False,excluded_resnames=excluded_resnames)
    single_altloc_atoms = ordered_atom_lookup.select_atoms_by(exclude_atom_names=["N","CA","C","O","H","H2","H3","HA"],
                                                altlocs=[parent_altlocs[0]])
    child_atom_tags = [DisorderedTag.from_atom(a) for a in single_altloc_atoms]

    split_model_path=os.path.join(out_dir,UntangleFunctions.model_handle(model_path)+"_split.pdb")
    split_specific(model_path,child_parent_altlocs_dict,child_atom_tags,out_path=split_model_path)

    text=create_all_child_restraints(model_path,child_parent_altlocs_dict,child_atom_tags)
    out_path=os.path.join(out_dir,"split_conformations_restraints.eff")
    with open(out_path,"w") as f:
        f.write(text)

if __name__ == "__main__":
    model_path="/home/speno/Untangler/output/2conf_start.pdb"
    child_parent_altlocs_dict={"C":"A"}
    run(model_path,child_parent_altlocs_dict)

# %%
