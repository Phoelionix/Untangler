#%%
import sys,os,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from ConformationTree.split_conformation_geoms import create_all_child_restraints
from LinearOptimizer.Tag import *
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler
from LinearOptimizer.Input import LP_Input
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from Bio.PDB import PDBParser,Structure,PDBIO
from UntangleFunctions import parse_symmetries_from_pdb,UNTANGLER_WORKING_DIRECTORY,prepare_pdb
from ConformationTree.split_pdb import split_specific


def run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=False):
    out_dir = os.path.join(UntangleFunctions.UNTANGLER_WORKING_DIRECTORY,"ConformationTree","output")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    single_altloc=sorted(list(UntangleFunctions.get_altlocs_from_pdb(model_path)[0]))[0]

    excluded_resnames=["CYS","GLY","PRO"]
    ordered_atom_lookup = OrderedAtomLookup(model_path, waters=False,excluded_resnames=excluded_resnames)
    atoms = ordered_atom_lookup.select_atoms_by(exclude_atom_names=["N","CA","C","O","H","H2","H3","HA"])
    child_atom_tags = list(set([DisorderedTag.from_atom(a) for a in atoms]))

    split_model_path=os.path.join(out_dir,UntangleFunctions.model_handle(model_path)+"_split.pdb")
    split_specific(model_path,child_parent_altlocs_dict,child_atom_tags,out_path=split_model_path,preserve_parent_altlocs=preserve_parent_altlocs)

    def get_out_path(model_handle,out_tag):
        output_dir = os.path.join(UNTANGLER_WORKING_DIRECTORY,"output","")
        return f"{output_dir}{model_handle}_{out_tag}.pdb"

    fmted_model = get_out_path(UntangleFunctions.model_handle(model_path),"fmtd")
    prepare_pdb(model_path,fmted_model,
                ring_name_grouping=False)

    text=create_all_child_restraints(fmted_model,child_parent_altlocs_dict,child_atom_tags)
    out_path=os.path.join(out_dir,"split_conformations_restraints.eff")
    with open(out_path,"w") as f:
        f.write(text)

if __name__ == "__main__":
    #model_path="/home/speno/Untangler/output/2conf_start.pdb"
    #model_path="/home/speno/Untangler/output/4PSS_split_refined.pdb"
    model_path="/home/speno/Untangler/data/4PSS_split.pdb"
    #child_parent_altlocs_dict={"C":"A"}
    # child_parent_altlocs_dict={"C":"A","D":"A","E":"A","F":"A","c":"B","d":"B","e":"B","f":"B"}
    # run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=False)
    child_parent_altlocs_dict={"C":"A","D":"B"}
    run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=True)

# %%
