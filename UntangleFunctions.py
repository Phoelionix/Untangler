


#%%
from Bio.PDB.Structure import Structure#, parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
# from typing import List # if your version of python is older than 3.9 uncomment this and replace instances of "list[" with "List[" 
import os
import subprocess
import matplotlib
from time import sleep
import numpy as np
matplotlib.use('AGG')
import shutil

ATOMS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'
       +' Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe').split()
NUM_E = {}
i = 0
for symbol in ATOMS:
    #TODO ions
    i+=1
    NUM_E[symbol] = i

UNTANGLER_WORKING_DIRECTORY= os.path.join(os.path.abspath(os.getcwd()),"")

def model_handle(model_path):
    assert model_path[-4:]==".pdb"
    return os.path.basename(model_path)[:-4]

def pdb_data_dir():
    return os.path.join(UNTANGLER_WORKING_DIRECTORY,"data","")

def separated_conformer_pdb_dir():
    return os.path.join(UNTANGLER_WORKING_DIRECTORY,"output","conformer_subsets","")


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

def get_score(score_file,phenixgeometry_only=False,verbose=True):
    assert os.path.exists(score_file)
    with open(score_file,'r') as f:
        for line in f: 
            if line.strip().strip('\n') == "":
                continue
            if verbose:
                print(line.strip('\n'))
            combined_score,wE_score,Rwork, Rfree = [float(v) for v in line.split()[:4]]
            if verbose:
                if phenixgeometry_only:
                    print(f"Python read wE (quick): {wE_score}")
                else:
                    print(f"Python read wE: {wE_score}")
            if verbose:
                print()
            return combined_score,wE_score,Rwork, Rfree


def batch_create_score_files(pdb_file_path,log_out_folder_path):
    assert False, "Unimplemented"

def create_score_file(pdb_file_path,log_out_folder_path,ignore_H=False,turn_off_cdl=False,reflections_for_R:str=None,skip_fail=False,timeout_mins=5): 
    # model_and_reflections_for_R overrides the R and R free values in the pdb path.

    holton_folder_path = UNTANGLER_WORKING_DIRECTORY+"StructureGeneration/"

    #generate_holton_data_shell_file=self.holton_folder_path+'GenerateHoltonData.sh'
    assert pdb_file_path[-4:]==".pdb"
    handle = os.path.basename(pdb_file_path)[:-4]
    #generate_holton_data_shell_file=self.holton_folder_path+'GenerateHoltonDataOriginal.sh' # for testing...
    generate_holton_data_shell_file=holton_folder_path+'GenerateHoltonData.sh'

    score_file=score_file_name(pdb_file_path,ignore_H=ignore_H)
    if os.path.exists(score_file):
        os.remove(score_file)

    assert os.path.exists(pdb_file_path)
    max_attempts=5
    i=0
    # Sometimes get seg faults
    while i < max_attempts and not os.path.exists(score_file):
        if i >0:
            print("Error, retrying")
        i+=1 
        sleep(1)
        print("--==--")
        rel_path = os.path.relpath(pdb_file_path,start=holton_folder_path)
        #with open(log_out_folder_path+f"{handle}_log.txt","w") as log:
        args = ["bash", f"{generate_holton_data_shell_file}",f"{rel_path}"]
        if ignore_H:
            args.append("-H")
        if turn_off_cdl:
            args.append("-C")
        if reflections_for_R is not None:
            Rwork,Rfree=get_R(pdb_file_path,reflections_for_R)
            args+= ["-W",str(Rwork),"F",str(Rfree)]
        print (f"|+ Running: {' '.join(args)}")
        try:
            subprocess.run(args,timeout=60*timeout_mins)#,stdout=log)
        except Exception as e:
            print(f"Error: {e}")
    if os.path.exists(score_file):
        print("finished")
        print("--==--")
    elif skip_fail:
        print("failed!")
    else:
        raise Exception(f"Failed to locate expected output file {score_file}")    


    return score_file

def copy_model_and_geo(source,dest):
    shutil.copy(source,dest)
    old_handle,new_handle = [model_handle(f) for f in (source,dest)]
    geo_dir=f"{UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/"

    for suffix in [".geo","_clashes.txt","_score.txt"]:
        shutil.copy(f"{geo_dir}{old_handle}{suffix}",f"{geo_dir}{new_handle}{suffix}")

def clear_geo(excluded_suffixes=["_start.geo","_fmtd.geo"]):
    holton_folder_path = os.path.join(UNTANGLER_WORKING_DIRECTORY,"StructureGeneration","")
    geo_path = os.path.join(holton_folder_path,"HoltonOutputs")
    for filename in os.listdir(geo_path):
        if filename.endswith(".geo"):
            for suffix in excluded_suffixes:
                 if filename.endswith(suffix):
                     break
            else: 
                os.remove(os.path.join(geo_path,filename))

def score_file_name(pdb_file_path,ignore_H=False):
    holton_folder_path = os.path.join(UNTANGLER_WORKING_DIRECTORY,"StructureGeneration","")
    assert pdb_file_path[-4:]==".pdb"
    handle = os.path.basename(pdb_file_path)[:-4]+("_ignoreH" if ignore_H else "")
    return holton_folder_path+f'HoltonOutputs/{handle}_score.txt'
def geo_file_name(pdb_file_path):
    holton_folder_path = os.path.join(UNTANGLER_WORKING_DIRECTORY,"StructureGeneration","")
    assert pdb_file_path[-4:]==".pdb"
    handle = os.path.basename(pdb_file_path)[:-4]
    return holton_folder_path+f'HoltonOutputs/{handle}.geo'

def assess_geometry_wE(pdb_file_path,log_out_folder_path,ignore_H=False,turn_off_cdl=False):
    score_file = create_score_file(pdb_file_path,log_out_folder_path,ignore_H=ignore_H,turn_off_cdl=turn_off_cdl)
    assert score_file == score_file_name(pdb_file_path,ignore_H=ignore_H)
    return get_score(score_file)

            
def res_is_water(res):
    return res.get_id()[0]!= " " or res.get_resname()=="HOH"

            

def get_R(pdb_file_path,reflections_path):
   

    print("--==--")
    #rel_path = os.path.relpath(pdb_file_path,start=holton_script_path)
    get_R_data_shell_file = f"{UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/GetRData.sh"
    args=["bash", f"{get_R_data_shell_file}",f"{pdb_file_path}",f"{reflections_path}"]
    print (f"|+ Running: {' '.join(args)}")
    subprocess.run(args)#,stdout=log)
    print("finished")
    print("--==--")

    log_file = f"{UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/output/SF_check_{model_handle(pdb_file_path)}.log"

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

def parse_symmetries_from_pdb(pdb_file_path):
    '''
    Parses the symmetry transformations from the pdb file into ndarrays and adds each to the crystal.  
    No unit operations are performed.
    Also returns the parsed matrices. 
    '''
    at_SYMOP = False
    at_symmetry_xformations = False
    sym_factor = np.zeros((3,3))
    sym_trans = np.zeros((3))
    sym_labels = []
    sym_mtces_parsed = []        
    with open(pdb_file_path) as pdb_file:
        for line in pdb_file:
            line = line.strip()
            # End data - Check if left section of data; flagged by the line containing solely "REMARK ###". 
            if len(line) <= 10: 
                if at_SYMOP:
                    at_SYMOP = False
                if at_symmetry_xformations: 
                    at_symmetry_xformations = False
            if line[0:4] == "ATOM":
                break
            # Symmetry operators (purely for terminal output - doesn't affect simulation)
            if at_SYMOP:
                sym_labels.append(line[17:21]+": ("+line[21:].strip()+")")
            # Add each symmetry to target
            if at_symmetry_xformations:
                entries = line.split()[2:]
                x = int(entries[0][-1]) - 1  
                entries = [float(f) for f in entries[2:]]
                # x/y index rows/cols
                for y, element in enumerate(entries[:-1]):
                    sym_factor[x,y] = element
                sym_trans[x] = entries[-1]
                if x == 2:
                    sym_mtces_parsed.append([sym_factor.copy(),sym_trans.copy()])


            # symmetry matrices
            if line == "REMARK 290 RELATED MOLECULES.":
                at_symmetry_xformations = True
    assert len(sym_mtces_parsed)!=0, f"Remark 290 symmetries appear to be missing from {pdb_file_path}"
    return sym_mtces_parsed

def get_sym_xfmed_point(R,symmetry):
    sym_rot, sym_trans = symmetry
    return (sym_rot @ R.T).T+ sym_trans

def two_key_read(dict_obj,dict_name,key1,key2):
    try: 
        return dict_obj[key1][key2]
    except:
        if key1 not in dict_obj:
            exception = f"{key1} not in {dict_name}"
            for k in dict_obj:
                print(k)
        else:
            exception = f"{key2} not in {dict_name}[{key1}]"
            for k in dict_obj[key1]:
                print(k)
        raise Exception(exception)
    

def prepare_pdb(pdb_path,out_path,sep_chain_format=False,altloc_from_chain_fix=False):
        # Gets into format we expect. !!!!!!Assumes single chain!!!!!
        def replace_occupancy(line,occ):
            occ=f"{occ:.3f}"
            occ = ' '*(6-len(occ))+occ
            return line[:54] + occ + line[60:]
        def replace_chain(line,chain_id):
            chain_id = str(chain_id)
            assert len(chain_id)==1
            return line[:21]+chain_id+line[22:]
        def replace_serial_num(line,serial_num):
            serial_num = str(serial_num)
            serial_num = ' '*(5-len(serial_num))+serial_num
            return line[:6]+serial_num+line[11:]
        def replace_altloc(line,altloc):
            return line[:16]+altloc+line[17:]
        def replace_res_num(line,res_num):
            res_num = str(res_num)
            res_num = ' '*(4-len(res_num))+res_num
            return line[:22]+res_num+line[26:]
            
        protein_altlocs = []
        solvent_altlocs = []
        with open(pdb_path) as I:
            max_resnum=0
            start_lines = []
            solvent_lines=[]
            end_lines = []
            atom_dict:dict[str,dict[str,dict[str,str]]] = {}  
            last_chain=None
            solvent_res_names=["HOH"]
            solvent_chain_id = "z"
            warned_collapse=False

            for line in I:
                if line.startswith("TER") or line.startswith("ANISOU"):
                    continue
                start_strs_considered = ["ATOM","HETATM"]
                for s in start_strs_considered:
                    if line.startswith(s):
                        break
                else: # Not modifying
                    if len(atom_dict)==0:
                        start_lines+=line
                    else:
                        end_lines += line
                    continue

                # Modifying
                name = line[12:16].strip()
                altloc = line[16]
                resname = line[17:20]
                space = line[20]
                chain = line[21]
                resnum = int(line[22:26])
                if altloc == ' ' and altloc_from_chain_fix:
                    altloc = chain
                    line = replace_altloc(line,altloc)
                if resname in solvent_res_names:
                    solvent_lines.append(replace_chain(line,solvent_chain_id))
                    if altloc not in solvent_altlocs:
                        solvent_altlocs.append(altloc) 
                    continue
                assert len(end_lines)==0
                
                # Non-solvent atoms
                if not sep_chain_format and not warned_collapse and chain != last_chain and last_chain is not None:
                    print("Warning: Multiple chains detected. Collapsing chains into single chain")
                    warned_collapse=True
                if resnum not in atom_dict:
                    atom_dict[resnum] = {}
                
                assert (altloc != ' '), line 

                if altloc not in atom_dict[resnum]:
                    atom_dict[resnum][altloc] = {}
                    
                    if altloc not in protein_altlocs:
                        protein_altlocs.append(altloc) 
                
                atom_dict[resnum][altloc][name]=line  
                max_resnum=max(resnum,max_resnum)
                last_chain = chain
                continue
                    
        n=0
        # Add non-solvent atoms
        if not sep_chain_format: # format for untangler stuff
            protein_chain_id = "A"
            for res_atom_dict in atom_dict.values():
                for altloc_atom_dict in res_atom_dict.values():
                    for line in altloc_atom_dict.values():
                        n+=1
                        modified_line = line
                        
                        modified_line = replace_occupancy(modified_line,
                            1/len(protein_altlocs)) # Set occupancies to all be same
                        modified_line = replace_chain(modified_line,protein_chain_id)
                        modified_line = replace_serial_num(modified_line,n)
                        start_lines.append(modified_line)
        else: # Note that lines for each chain need to be contiguous in the file
            chain_dict={}
            for res_atom_dict in atom_dict.values():
                for altloc, altloc_atom_dict in res_atom_dict.items():
                    protein_chain_id = altloc
                    for line in altloc_atom_dict.values():
                        n+=1
                        modified_line = line
                        modified_line = replace_occupancy(modified_line,
                            1/len(protein_altlocs)) # Set occupancies to all be same
                        modified_line = replace_chain(modified_line,protein_chain_id)
                        modified_line = replace_serial_num(modified_line,n)
                        if altloc not in chain_dict:
                            chain_dict[altloc]=[]
                        chain_dict[altloc].append(modified_line)
            for _, lines in chain_dict.items():
                for modified_line in lines:
                    start_lines.append(modified_line)

        # Make sure waters don't share residue numbers with protein
        min_solvent_resnum=99999999
        for line in solvent_lines:
            solvent_resnum=int(line[22:26])
            min_solvent_resnum = min(solvent_resnum,min_solvent_resnum)
        shift = max_resnum-min_solvent_resnum + 1
        new_solvent_resnum_dict = {}
        for line in solvent_lines:
            solvent_resnum=int(line[22:26])
            # In case of gaps...
            if solvent_resnum not in new_solvent_resnum_dict:
                if (solvent_resnum+shift)-max_resnum > 1:
                    shift = max_resnum-solvent_resnum + 1 
                elif (solvent_resnum+shift)-max_resnum < 0:
                    print(f"Warning, are solvent res nums out of order? {max_resnum,solvent_resnum,shift, min_solvent_resnum, max_resnum, out_path}")
                max_resnum = shift+solvent_resnum
                new_solvent_resnum_dict[solvent_resnum]=max_resnum
            
            
            modified_line = replace_res_num(line,new_solvent_resnum_dict[solvent_resnum])
            start_lines.append(modified_line)

        with open(out_path,'w') as O:
            O.writelines(start_lines+end_lines)


def get_altlocs_from_pdb(pdb_path):
            
        protein_altlocs = []
        solvent_altlocs = []
        with open(pdb_path) as I:
            solvent_res_names=["HOH"]

            for line in I:
                if line.startswith("TER") or line.startswith("ANISOU"):
                    continue
                start_strs_considered = ["ATOM","HETATM"]
                for s in start_strs_considered:
                    if line.startswith(s):
                        break
                else: 
                    continue 

                # Atom entries
                altloc = line[16]
                resname = line[17:20]
                if resname in solvent_res_names:
                    if altloc not in solvent_altlocs:
                        solvent_altlocs.append(altloc) 
                    continue
                
                # Non-solvent atoms
                
                assert (altloc != ' '), line 


                if altloc not in protein_altlocs:
                    protein_altlocs.append(altloc) 

        return set(protein_altlocs),set(solvent_altlocs)
# %%
