#%%
import sys,os,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from ConformationTree.split_conformation_geoms import create_all_child_restraints
from LinearOptimizer.Tag import *
from LinearOptimizer.RestraintsHandler import RestraintsHandler
from LinearOptimizer.Input import LP_Input
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from Bio.PDB import PDBParser,Structure,PDBIO
from UntangleFunctions import parse_symmetries_from_pdb, UNTANGLER_WORKING_DIRECTORY, prepare_pdb
from ConformationTree.split_pdb import split_specific


def run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=False,nonexistent_parent_from_child_priority_dict={},equalize_output_occupancies=False,never_alter_lone_waters=False, nonexistent_parents_replace_child=True,
        shake_new_conformers=0.1 # Angstrom
        ):
    out_dir = os.path.join(UntangleFunctions.UNTANGLER_WORKING_DIRECTORY,"ConformationTree","output")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    def get_out_path(model_handle,out_tag):
        output_dir = os.path.join(UNTANGLER_WORKING_DIRECTORY,"output","")
        return f"{output_dir}{model_handle}_{out_tag}.pdb"

    fmted_model = get_out_path(UntangleFunctions.model_handle(model_path),"fmtd")
    prepare_pdb(model_path,fmted_model,
            ring_name_grouping=False,
            allow_no_altloc=True)

    #excluded_resnames=["CYS","GLY","PRO"]
    excluded_resnames=[]
    ordered_atom_lookup = OrderedAtomLookup(fmted_model, waters=True,excluded_resnames=excluded_resnames)
    #atoms_to_split = ordered_atom_lookup.select_atoms_by(altlocs=[' '])
    atoms_to_split = ordered_atom_lookup.ordered_atoms
    site_tags = [DisorderedTag.from_atom(a) for a in atoms_to_split]
    assert len(site_tags)>0

    split_model_path=os.path.join(out_dir,UntangleFunctions.model_handle(model_path)+"_split_unrefined.pdb")
    split_specific(fmted_model,child_parent_altlocs_dict,site_tags,out_path=split_model_path,preserve_parent_altlocs=preserve_parent_altlocs,
                   split_waters=True,never_alter_lone_waters=never_alter_lone_waters,nonexistent_parent_from_child_priority_dict=nonexistent_parent_from_child_priority_dict,
                   nonexistent_parents_replace_child=nonexistent_parents_replace_child,shake_new=shake_new_conformers)
    
    if equalize_output_occupancies:
        prepare_pdb(split_model_path,split_model_path,
                    even_split_protein_occupancies=True,treat_solvent_identically_to_protein=True)
    print(f"Successfully written to {split_model_path}")

if __name__ == "__main__":
    model_path="/home/speno/Untangler/data/lys/5KXL.pdb"
    equalize_output_occupancies=True
    preserve_parent_altlocs=False
    nonexistent_parents_replace_child=True

    SING=' '
    #child_parent_altlocs_dict={}; nonexistent_parent_from_child_priority_dict={"A":"CDB"+SING,"B":"DCA"+SING,"C"+SING:"ABD ","D":"BAC"+SING}; nonexistent_parents_replace_child=False # Make contiguous conformations from 4-conf qFit. 
    child_parent_altlocs_dict={};nonexistent_parent_from_child_priority_dict={"A":"BC"+SING,"B":"CA"+SING,"C":"AB"+SING,"D":'BAC'+SING,"E":'CAB'+SING}; nonexistent_parents_replace_child=False # Make 5 contiguous conformations from 3-conf qFit. Could be more thought out.


    #child_parent_altlocs_dict={"A":" ","B":" ", "C":" ", "D":" ", "E": " "}; preserve_parent_altlocs=False

    #nonexistent_parent_from_child_priority_dict={" ":"ABCDE"}
    #child_parent_altlocs_dict={"C":"A","D":"A","E":"A","F":"A","c":"B","d":"B","e":"B","f":"B"}
    #child_parent_altlocs_dict={"A":" ","B":" ","C":" ","D":" ","E":" ","F":" "}
    #child_parent_altlocs_dict={"D":"A","E":"B","F": "C"}
    shake_new_conformers=0.5 # angstrom
    run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=preserve_parent_altlocs,equalize_output_occupancies=equalize_output_occupancies,
        nonexistent_parents_replace_child=nonexistent_parents_replace_child,
        nonexistent_parent_from_child_priority_dict=nonexistent_parent_from_child_priority_dict,
        shake_new_conformers=shake_new_conformers)

# rename to e.g. XXXX_Nconf_unrefined.pdb
#  bash /home/speno/Untangler/Refinement/Refine.sh /home/speno/Untangler/data/XXXX_Nconf_unrefined.pdb  data/XXXX.mtz -c 1.0 -u 1.0 -n 6 -q 0.1 -N -s 0.5

# %%
