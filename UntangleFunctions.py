


#%%
from Bio.PDB.Structure import Structure#, parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
# from typing import List # if your version of python is older than 3.9 uncomment this and replace instances of "list[" with "List[" 
import os
import subprocess
import matplotlib
matplotlib.use('AGG')



UNTANGLER_WORKING_DIRECTORY= os.path.join(os.path.abspath(os.getcwd()),"")

def pdb_data_dir():
    return UNTANGLER_WORKING_DIRECTORY+"data/"

def get_header(pdb_file_path):
    header = ""
    with open(pdb_file_path,'r') as f:  # TODO try just saving header to structure
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                break
            header += line
    return header

def save_structure(structure:Structure,header_reference_file_path,out_path,comments_to_write=None):
    io = PDBIO()

    io.set_structure(structure)
    out_folder_path = os.path.dirname(out_path)
    tmp_path = out_folder_path+"structure_combined_tmp.pdb"
    io.save(tmp_path)

    header = ""
    if comments_to_write is not None:
        header=comments_to_write
    header += get_header(header_reference_file_path)

    with open(tmp_path,'r') as f:
        with open(out_path,'w') as f2: 
            f2.write(header)
            f2.write(f.read())
    os.remove(tmp_path)

def get_score(score_file,phenixgeometry_only=False):
    with open(score_file,'r') as f:
        for line in f: 
            if line.strip().strip('\n') != "":
                print(line.strip('\n'))
                combined_score,wE_score,Rwork, Rfree = [float(v) for v in line.split()[:4]]
                if phenixgeometry_only:
                    print(f"Python read wE (quick): {wE_score}")
                else:
                    print(f"Python read wE: {wE_score}")
                print()
                return combined_score,wE_score,Rwork, Rfree


def create_score_file(log_out_folder_path,pdb_file_path,phenixgeometry_only=False):
    holton_folder_path = UNTANGLER_WORKING_DIRECTORY+"StructureGeneration/"

    #generate_holton_data_shell_file=self.holton_folder_path+'GenerateHoltonData.sh'
    assert pdb_file_path[-4:]==".pdb"
    handle = os.path.basename(pdb_file_path)[:-4]
    #generate_holton_data_shell_file=self.holton_folder_path+'GenerateHoltonDataOriginal.sh' # for testing...
    generate_holton_data_shell_file=holton_folder_path+'GenerateHoltonData.sh'
    if phenixgeometry_only:
        generate_holton_data_shell_file=holton_folder_path+'GenerateHoltonDataQuick.sh' #

    print("--==--")
    print (f"running {generate_holton_data_shell_file}")
    rel_path = os.path.relpath(pdb_file_path,start=holton_folder_path)
    with open(log_out_folder_path+f"{handle}_log.txt","w") as log:
        args = ["bash", f"{generate_holton_data_shell_file}",f"{rel_path}"]
        print(f"Running {args}")
        subprocess.run(args,stdout=log)
    print("finished")
    print("--==--")
    return score_file_name(pdb_file_path)

def score_file_name(pdb_file_path):
    holton_folder_path = UNTANGLER_WORKING_DIRECTORY+"StructureGeneration/"
    assert pdb_file_path[-4:]==".pdb"
    handle = os.path.basename(pdb_file_path)[:-4]
    return holton_folder_path+f'HoltonOutputs/{handle}_score.txt'

def assess_geometry_wE(pdb_file_path,log_out_folder_path,phenixgeometry_only=False):
    score_file = create_score_file(log_out_folder_path,pdb_file_path,phenixgeometry_only)
    assert score_file == score_file_name(pdb_file_path)
    return get_score(score_file,phenixgeometry_only=phenixgeometry_only)

            
def res_is_water(res):
    return res.get_id()[0]!= " " or res.get_resname()=="HOH"

            

def get_R(pdb_file_path,reflections_path):
    assert pdb_file_path[-4:]==".pdb"
    model_handle = os.path.basename(pdb_file_path)[:-4]

    print("--==--")
    #rel_path = os.path.relpath(pdb_file_path,start=holton_script_path)
    get_R_data_shell_file = f"{UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/GetRData.sh"
    args=["bash", f"{get_R_data_shell_file}",f"{pdb_file_path}",f"{reflections_path}"]
    print (f"$Running {args}")
    subprocess.run(args)#,stdout=log)
    print("finished")
    print("--==--")

    log_file = f"{UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/output/SF_check_{model_handle}.log"

    Rwork=Rfree=None
    with open(log_file,'r') as f:
        for line in f:
            if line.startswith("  r_work:"):
                Rwork = float(line.split(" ")[-1])
            if line.startswith("  r_free:"):
                Rfree= float(line.split(" ")[-1])
                assert Rwork!=None
                break
    print(Rwork,Rfree)
    return Rwork,Rfree

# nonH_fullname_list is with whitespaces
# quite possibly the worst parser ever written.
def H_get_parent_fullname(H_name:str,nonH_fullname_list:list[str],debug_print=False): 
        parent_name=None
        H_name = H_name.strip()
        basic_dict = dict(
            H = " N  ",
            H1 = " N  ",
            H2 = " N  ",
            H3 = " N  ",
            HA = " CA ",
            HB = " CB ",
        )
        # determine parent name
        if H_name in basic_dict:
            parent_name = basic_dict[H_name]
        else:
            identifier = H_name[1:]
            for k in nonH_fullname_list:
                if k.strip()[0] == "H":
                   raise Exception("nonH atom list contains an atom name starting with H!") 
                strip_digit=False
                if identifier[-1].isnumeric(): # This is very annoying
                    if not k.strip()[-1].isnumeric():
                        strip_digit=identifier[-1].isnumeric()
                    else:
                        strip_digit=identifier[-2].isnumeric()
                condition  = k[2:].strip()==identifier
                if strip_digit:
                    condition  = k[2:].strip()==identifier[:-1]
                if condition:
                    #print(H.name, identifier)
                    if parent_name is not None:
                        #if debug_print:
                        print(H_name, identifier)
                        for k in nonH_fullname_list:
                            print(k,k[2:].strip())
                        assert False
                    parent_name = k
        if parent_name is None:
            if debug_print:
                print(H_name, identifier)
                for k in nonH_fullname_list:
                    print(k,k[2:].strip())
                assert False
        return parent_name

# %%
