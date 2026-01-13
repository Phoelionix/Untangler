import sys,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from LinearOptimizer.Tag import *
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler
from LinearOptimizer.Input import LP_Input
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from Bio.PDB import PDBParser,Structure,PDBIO
from UntangleFunctions import parse_symmetries_from_pdb





def split_specific(pdb_path,child_parent_altlocs_dict,child_atom_tags:list[DisorderedTag],out_path=None,
                   sep_chain_format=False,protein_altloc_from_chain_fix=False,missing_water_altloc_fix=True):
    # Splits conformers of atoms (atoms_to_split) according to child_parent_altlocs_dict


    assert pdb_path[-4:] == ".pdb"
    if out_path is None:
        out_path = pdb_path[:-4]+"_split.pdb"
    
    
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
    def replace_occupancy(line,occ):
        occ=f"{occ:.3f}"
        occ = ' '*(6-len(occ))+occ
        return line[:54] + occ + line[60:]
        
    protein_altlocs = []
    #solvent_altlocs = []
    with open(pdb_path) as I:
        max_resnum=0
        max_serial_num=0
        start_lines = []
        solvent_lines=[]
        end_lines = []
        atom_dict:dict[str,dict[str,dict[str,str]]] = {}  
        last_chain=None
        solvent_res_names=["HOH"]
        solvent_chain_id = "z"
        warned_collapse=False

        convert_to_ensemble = True
        lines = I.readlines()
        for line in lines: 
            if line.startswith("TER") or line.startswith("ANISOU"):
                continue
            start_strs_considered = ["ATOM","HETATM"]
            for s in start_strs_considered:
                if line.startswith(s):
                    altloc = line[16]
                    if altloc!=' ':
                        convert_to_ensemble=False
        
        resnum_chain_append_mod=0 
        last_resnum=None
        for line in lines:
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
            if resname in solvent_res_names:
                solvent_lines.append(line)  # Modified further below
                continue

            if altloc == ' ':
                if convert_to_ensemble:
                    altloc = "A" 
                elif protein_altloc_from_chain_fix:
                    altloc = chain
                if altloc != ' ' and altloc not in protein_altlocs:
                    protein_altlocs.append(altloc)
                line = replace_altloc(line,altloc)
            assert len(end_lines)==0
                
            if not sep_chain_format and chain != last_chain and last_chain is not None:
                if not warned_collapse:
                    print("Warning: Multiple chains detected. Collapsing chains into single chain")
                if altloc != chain: # XXX
                    if not warned_collapse:
                        print("Appending resnums")
                    this_line_resnum=int(line[22:26])
                    resnum_chain_append_mod=last_resnum+1-this_line_resnum
                warned_collapse=True

            last_resnum=resnum = int(line[22:26])+resnum_chain_append_mod
            line = replace_res_num(line,resnum)

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

    #assert len(set(protein_altlocs) & set(solvent_altlocs))==len(protein_altlocs)==len(solvent_altlocs), (protein_altlocs,solvent_altlocs)

    chain_dict={}
    for res_atom_dict in atom_dict.values():
        #for d in range(2): # Split loop
        for parent_altloc, altloc_atom_dict in res_atom_dict.items(): 
            child_altlocs=[k for k,v in child_parent_altlocs_dict.items() if parent_altloc in v]
            num_relevant_altlocs=1+len(child_altlocs)
            for d in range(num_relevant_altlocs):
                for line in altloc_atom_dict.values():
                    if sep_chain_format:
                        protein_chain_id=altloc
                    else:
                        protein_chain_id="A"
                    parent_altloc=line[16]
                    max_serial_num+=1
                    modified_line = replace_serial_num(line,max_serial_num)
                    modified_line = replace_chain(modified_line,protein_chain_id)
                    atom_name=line[12:16].strip()
                    res_num = int(line[22:26])
                    site_tag = DisorderedTag(res_num,atom_name)
                    atom_being_split = site_tag in child_atom_tags
                    if atom_being_split:
                        occupancy=float(line[54:60])
                        modified_line = replace_occupancy(modified_line,occupancy/num_relevant_altlocs)
                        if d > 0:
                            modified_line = replace_altloc(modified_line,child_altlocs[d-1])
                    elif d>0:
                        continue
                    if protein_chain_id not in chain_dict:
                        chain_dict[protein_chain_id]=[]
                    chain_dict[protein_chain_id].append(modified_line)

    for _, lines in chain_dict.items():
        for modified_line in lines:
            start_lines.append(modified_line)

    # Make sure waters don't share residue numbers with protein
    min_solvent_resnum=99999999
    for line in solvent_lines:
        solvent_resnum=int(line[22:26])
        min_solvent_resnum = min(solvent_resnum,min_solvent_resnum)
        shift = max_resnum-min_solvent_resnum + 1
    for line in solvent_lines:
        solvent_resnum=int(line[22:26])
        modified_line = line
        modified_line = replace_res_num(modified_line,shift+solvent_resnum)
        modified_line=replace_serial_num(modified_line,max_serial_num)
        max_serial_num+=1
        split_waters=False
        if split_waters:
            # start_lines.append(modified_line)
            # splt = replace_altloc(modified_line,split_altloc(modified_line))
            # splt = replace_serial_num(splt,max_serial_num)
            # max_serial_num+=1
            # start_lines.append(splt)
            assert False
        else:
            start_lines.append(modified_line)


    with open(out_path,'w') as O:
        O.writelines(start_lines+end_lines)



