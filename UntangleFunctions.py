


#%%
from Bio.PDB.Structure import Structure#, parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom
# from typing import List # if your version of python is older than 3.9 uncomment this and replace instances of "list[" with "List[" 
import os
import subprocess
from time import sleep
import numpy as np
import shutil
import itertools

ATOMS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'
       +' Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe').split()
NUM_E = {}
i = 0
for symbol in ATOMS:
    #TODO ions
    i+=1
    NUM_E[symbol] = i

UNTANGLER_WORKING_DIRECTORY= os.path.join(os.path.abspath(os.getcwd()),"")

def model_handle(model_handle_or_path):
    if model_handle_or_path[-4:]==".pdb":
        return os.path.basename(model_handle_or_path)[:-4]
    else: 
        # assume is already a handle
        assert os.path.sep not in model_handle_or_path
        return model_handle_or_path


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



def create_score_file(pdb_file_path,ignore_H=False,turn_off_cdl=False,reflections_for_R:str=None,skip_fail=False,timeout_mins=5): 
    # model_and_reflections_for_R overrides the R and R free values in the pdb path.

    holton_folder_path = UNTANGLER_WORKING_DIRECTORY+"StructureGeneration/"

    #generate_holton_data_shell_file=self.holton_folder_path+'GenerateHoltonData.sh'
    assert pdb_file_path[-4:]==".pdb"
    handle = os.path.basename(pdb_file_path)[:-4]
    #generate_holton_data_shell_file=self.holton_folder_path+'GenerateHoltonDataOriginal.sh' # for testing...
    generate_holton_data_shell_file=holton_folder_path+'GenerateHoltonData.sh'

    score_file=score_file_name(pdb_file_path,ignore_H=ignore_H,turn_off_cdl=turn_off_cdl)
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

def score_file_name(model_handle_or_path,ignore_H=False,turn_off_cdl=False):
    handle = model_handle(model_handle_or_path)+("_ignoreH" if ignore_H else "")+("_noCDL" if turn_off_cdl else "")
    return os.path.join(UNTANGLER_WORKING_DIRECTORY,"StructureGeneration",'HoltonOutputs',f'{handle}_score.txt')

def geo_file_name(model_handle_or_path):
    handle = model_handle(model_handle_or_path)
    return os.path.join(UNTANGLER_WORKING_DIRECTORY,"StructureGeneration","HoltonOutputs",f"{handle}.geo")

def assess_geometry_wE(pdb_file_path,ignore_H=False,turn_off_cdl=False):
    score_file = create_score_file(pdb_file_path,ignore_H=ignore_H,turn_off_cdl=turn_off_cdl)
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
def H_get_parent_fullname(H_name:str,nonH_fullname_list:list[str],debug_print=True): 
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
    


def relabel_ring(pdb_path):

    def sqrdist(a:Atom,b:Atom):
            #print(a_shifted,b_shifted)
            return np.sum((a.get_coord()-b.get_coord())**2)

    # TODO do with NQH?

    struct = PDBParser().get_structure("struct",pdb_path)
    # Finds labels such that RMSD for the two groups is as low as possible
    indistinguish_dict = {"CD1":"CD", 
                          "CD2":"CD",
                          "CE1":"CE",
                          "CE2":"CE"}
    
    anchors={"CD1":"CE1", "CD2":"CE2"}
    riding_hydrogens={"CD1":"HD1", "CD2":"HD2","CE1":"HE1", "CE2":"HE2"}
    original_label_dict={}
    ring_relabel_dict:dict[int,str]={}
    total_sd_original=0
    total_sd_after=0

    def angle(a,b,c):
        v1 = a - b
        v2 = c - b
        def unit_vector(vector):
            return vector / np.linalg.norm(vector)
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        return np.arccos(np.dot(v1_u, v2_u))*180/np.pi

    for residue in [res for res in struct.get_residues() if res.get_resname() in ["TYR","PHE"]]:
        ring_group_dict = {key:[] for key in ["CD","CE"]}
        avg_CG_coord=None
        #CG_coords={}
        hydrogen_sn_dict={}
        for Datom in residue.get_atoms():
            for atom in Datom:
                name = atom.get_name()
                if name in indistinguish_dict:
                    ring_group_dict[indistinguish_dict[name]].append(atom)
                    original_label_dict[atom.get_serial_number()]=name
                elif name in riding_hydrogens.values():
                    if name not in hydrogen_sn_dict:
                        hydrogen_sn_dict[name]={}
                    hydrogen_sn_dict[name][atom.get_altloc()]=atom.get_serial_number()
            if Datom.get_name()=="CG":
                avg_CG_coord=np.mean([a.get_coord() for a in Datom],axis=0)
                #CG_coords[atom.get_altloc()]=atom.get_coord()
        anchor_vals={"CD1":None,"CD2":None}; groups=["CD","CE"]
        for group_kind in groups:
            atoms=ring_group_dict[group_kind]
            group_names= [key for key in indistinguish_dict if indistinguish_dict[key]==group_kind]
            assert len(group_names)==2

            #grouping_options:itertools.combinations[tuple[list[Atom],list[Atom]]] = itertools.combinations(atoms, int(len(atoms)/2))


            # for all 2 element combinations of set
            grouping_options=[]
            for perm in itertools.combinations(atoms,int(len(atoms)/2)):
                pair = (perm,tuple([a for a in atoms if a not in perm]))
                # need same number of altlocs in each.
                altlocs_a,altlocs_b = [a.get_altloc() for a in pair[0]], [b.get_altloc() for b in pair[1]]
                if set(altlocs_a)!=set(altlocs_b):
                    continue
                pair = set(pair)
                if pair not in grouping_options:
                    grouping_options.append(pair)

            #print(atoms)
            assert len(grouping_options)>0

    
            best_option:tuple[list[Atom],list[Atom]] = None
            lowest_sd=np.inf
            best_means=None
            for a,b in grouping_options:
                deviations=[]
                means=[]
                for g in a,b:
                    mean = np.mean([atom.get_coord() for atom in g],axis=0)
                    # if residue.get_id()[1]==47:
                    #     print([(atom.get_name(),atom.get_coord()) for atom in g] )
                    #     print(mean)
                    deviations.extend([atom.get_coord()-mean for atom in g])
                    means.append(mean)
                rmsd = np.sqrt(np.mean([np.sum(dev**2) for dev in deviations]))
                # if residue.get_id()[1]==47:
                #     print(rmsd)
                if rmsd < lowest_sd:
                    lowest_sd=rmsd
                    best_option=(a,b)
                    best_means=means
                if len(set([atom.get_name() for atom in a]))==1:
                    assert len(set([atom.get_name() for atom in b]))==1
                    original_sd=rmsd
            assert best_option is not None, grouping_options

            need_anchor=False
            for group_name, mean in zip(group_names,best_means):
                #print(group_name)
                if group_name in anchor_vals:
                    anchor_vals[group_name]=mean
                else:
                    need_anchor=True
                    assert np.all([v is not None for v in anchor_vals.values()])
                #print(anchor_vals)
            anchor_names = list(anchor_vals)
            final_names_dict={group_names[i]: best_option[i] for i in range(2) }
            if need_anchor:
                order_a = (anchor_names,{anchor_names[0]:best_means[0],anchor_names[1]:(best_means[1])})
                order_b = (list(reversed(anchor_names)),{anchor_names[1]:best_means[0],anchor_names[0]:(best_means[1])})
                
                best=np.inf
                angles=[]
                ideal_angle=180
                for i, (mod_anchor_names, mean_dict) in enumerate((order_a,order_b)):
                    for j, (anchor_name, mean) in enumerate(mean_dict.items()):
                        angles.append(angle(mean,anchor_vals[anchor_name],avg_CG_coord))
                        deviations.append(angles[-1]-ideal_angle)
                        # if residue.get_id()[1]==15:
                        #     print(np.sqrt(np.sum((mean-anchor_vals[anchor_name])**2)),np.sqrt(np.sum((avg_CG_coord-anchor_vals[anchor_name])**2)))
                        #     print(anchor_name,best_option[j],mean,angles[-1])
                        #     print(mean,anchor_vals[anchor_name],avg_CG_coord)
                    rmsd = np.sqrt(np.mean([np.sum(dev**2) for dev in deviations]))
                    if rmsd<best:
                        best=rmsd
                        final_names_dict = {anchors[n]:g for n,g in zip(mod_anchor_names,best_option)}
            # if residue.get_id()[1]==47:
            #     print(final_names_dict)
            # TODO CD2 and CE2 correspond to covalent bond. So should decide whether label is 1 or 2 
            
            
            assert len(final_names_dict)==len(best_option)
            for new_atom_name, atoms in final_names_dict.items():
                for atom in atoms:
                    assert atom.get_serial_number() not in ring_relabel_dict
                    ring_relabel_dict[atom.get_serial_number()]=new_atom_name
                    for hydrogen_name in hydrogen_sn_dict:
                        for hydrogen_altloc, hydrogen_sn in hydrogen_sn_dict[hydrogen_name].items():
                            if riding_hydrogens[atom.get_name()]==hydrogen_name and hydrogen_altloc==atom.get_altloc():
                                ring_relabel_dict[hydrogen_sn]=riding_hydrogens[new_atom_name]

            total_sd_original+=original_sd 
            total_sd_after+=lowest_sd

    #print(original_label_dict)
    #print(ring_relabel_dict)
    print(f"{total_sd_original:.2f}>>{total_sd_after:.2f}")
    return ring_relabel_dict
    

def prepare_pdb(pdb_path,out_path,sep_chain_format=False,altloc_from_chain_fix=False,ring_name_grouping=False):
        # Gets into format we expect. !!!!!!Assumes single chain!!!!!
        # Relabels ring atoms CE1/CE2, CD1/CD2 so that all with same label are closest         
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
        def replace_name(line,name):
            if len(name)==3:
                name = ' '+name
            elif len(name)==2:
                name = ' '+name+' '
            elif len(name)==1:
                name = ' '+name+'  '
            else:
                assert len(name)==4
            return line[:12]+name+line[16:]
            
        
        if ring_name_grouping:
            ring_relabel_dict = relabel_ring(pdb_path)

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

                if ring_name_grouping:
                    if resname in ["TYR","PHE"] and name in ["CD1","CD2","CE1","CE2","HD1","HD2","HE1","HE2"]:
                        original_serial_num=int(line[6:11])
                        name = ring_relabel_dict[original_serial_num]
                        line = replace_name(line,name)
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
