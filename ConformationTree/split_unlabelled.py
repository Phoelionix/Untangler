#%%
import sys,os,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from ConformationTree.split_conformation_geoms import create_all_child_restraints
from LinearOptimizer.Tag import *
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler
from LinearOptimizer.Input import LP_Input
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from Bio.PDB import PDBParser,Structure,PDBIO
from UntangleFunctions import parse_symmetries_from_pdb, UNTANGLER_WORKING_DIRECTORY, prepare_pdb
from ConformationTree.split_pdb import split_specific


def run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=False,nonexistent_parent_from_child_priority_dict={}):
    out_dir = os.path.join(UntangleFunctions.UNTANGLER_WORKING_DIRECTORY,"ConformationTree","output")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    def get_out_path(model_handle,out_tag):
        output_dir = os.path.join(UNTANGLER_WORKING_DIRECTORY,"output","")
        return f"{output_dir}{model_handle}_{out_tag}.pdb"

    fmted_model = get_out_path(UntangleFunctions.model_handle(model_path),"fmtd")
    prepare_pdb(model_path,fmted_model,
            ring_name_grouping=False)

    #excluded_resnames=["CYS","GLY","PRO"]
    excluded_resnames=[]
    ordered_atom_lookup = OrderedAtomLookup(fmted_model, waters=True,excluded_resnames=excluded_resnames)
    #atoms_to_split = ordered_atom_lookup.select_atoms_by(altlocs=[' '])
    atoms_to_split = ordered_atom_lookup.ordered_atoms
    site_tags = [DisorderedTag.from_atom(a) for a in atoms_to_split]
    assert len(site_tags)>0

    split_model_path=os.path.join(out_dir,UntangleFunctions.model_handle(fmted_model)+"_split.pdb")
    split_specific(model_path,child_parent_altlocs_dict,site_tags,out_path=split_model_path,preserve_parent_altlocs=preserve_parent_altlocs,
                   split_waters=True,nonexistent_parent_from_child_priority_dict=nonexistent_parent_from_child_priority_dict,
                   force_lone_altloc_label=" ")


if __name__ == "__main__":
    #model_path="/home/speno/Untangler/data/4PSS.pdb"
    model_path="/home/speno/Untangler/data/cov63/cov63_2confR.pdb"
    #child_parent_altlocs_dict={"C":"A"}
    # child_parent_altlocs_dict={"C":"A","D":"A","E":"A","F":"A","c":"B","d":"B","e":"B","f":"B"}
    # run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=False)
    child_parent_altlocs_dict={"A":" ","B":" "}
    nonexistent_parent_from_child_priority_dict={"A":"CD","B":"DC"}
    run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=False,nonexistent_parent_from_child_priority_dict=nonexistent_parent_from_child_priority_dict)

# %%
