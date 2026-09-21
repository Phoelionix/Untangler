#%%
import sys,os,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from ConformationTree.split_all import run

assert len(sys.argv)==3, "Usage: python3 split_qfit_to_6conf input_model.pdb output_model.pdb"
model_path,out_path=sys.argv[1:3]
equalize_output_occupancies=True
preserve_parent_altlocs=False
nonexistent_parents_replace_child=True
SING=' '
child_parent_altlocs_dict={};nonexistent_parent_from_child_priority_dict={"A":"BC"+SING,"B":"CA"+SING,"C":"AB"+SING,"D":'BAC'+SING,"E":'CAB'+SING,'F':'ACB'+SING}; nonexistent_parents_replace_child=False # Make 6 contiguous conformations from 3-conf qFit. Could be more thought out.

shake_new_conformers=0.5 # angstrom
run(model_path,child_parent_altlocs_dict,preserve_parent_altlocs=preserve_parent_altlocs,equalize_output_occupancies=equalize_output_occupancies,
    nonexistent_parents_replace_child=nonexistent_parents_replace_child,
    nonexistent_parent_from_child_priority_dict=nonexistent_parent_from_child_priority_dict,
    shake_new_conformers=shake_new_conformers,
    out_path=out_path)