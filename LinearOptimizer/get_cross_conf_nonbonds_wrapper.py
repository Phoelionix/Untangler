import subprocess
import os, sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from UntangleFunctions import UNTANGLER_WORKING_DIRECTORY
# Interface for phenix.python call
def get_cross_conf_nonbonds(pdb_file_path,use_previous=False,use_previous_nonbonds=False,nonbonds=False):
    # use_previous_nonbonds is useful, as the code gets all nonbond distances for atoms within a wide range of each other. 
    # But clashes need to be regenerated each time due to calling the probe program.

    kind='nonbonds' if nonbonds else 'clashes'
    tmp_file_path=f"tmp_{kind}_data.txt"
    if not use_previous and not (use_previous_nonbonds and nonbonds):
        if os.path.exists(tmp_file_path):
            os.remove(tmp_file_path)
        script_path=os.path.join(UNTANGLER_WORKING_DIRECTORY,"PhenixEnvScripts",f"cross_conformation_{kind}.py")
        args=["phenix.python",script_path,pdb_file_path,tmp_file_path]
        print (f"$Running {' '.join(args)}")
        subprocess.run(args)#,stdout=log)
    else:
        print(f"Warning: Using previous {kind} data files")


    out_items = []
    with open(tmp_file_path) as f:
        for line in f:
            pdb1,pdb2,vdw,is_symm = line.rstrip().split('|')
            vdw = float(vdw)
            is_symm = (is_symm.lower() == "true")
            out_items.append((pdb1,pdb2,vdw,is_symm))
    #os.remove(tmp_file_path)
    return out_items

if __name__ == "__main__":
    get_cross_conf_nonbonds(*sys.argv[1:])
