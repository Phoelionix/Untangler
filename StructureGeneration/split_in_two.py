
import sys, os





def split(pdb_path,out_path,sep_chain_format=False,protein_altloc_from_chain_fix=False,missing_water_altloc_fix=True):
    new_altloc_options = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    
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
    solvent_altlocs = []
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
            occupancy=float(line[54:60])
            line = replace_occupancy(line,occupancy/2)
            if resname in solvent_res_names:
                solvent_lines.append(line)  # Modified further below
                continue

            if altloc == ' ' and protein_altloc_from_chain_fix:
                altloc = chain
                line = replace_altloc(line,altloc)
            assert len(end_lines)==0
                
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
    new_solvent_lines=[]
    for line in solvent_lines:
        sublines= [line]
        altloc = line[16]
        if altloc == ' ' and missing_water_altloc_fix:
            sublines= []
            for a in protein_altlocs:
                occupancy=float(line[54:60])
                solv_line = replace_occupancy(line,occupancy/2)
                sublines.append(replace_altloc(solv_line,a))
        for solv_line in sublines:
            a = solv_line[16]
            #solv_line = replace_res_num(solv_line,max_resnum+1)
            #print(solv_line)
            #max_resnum+=1
            new_solvent_lines.append(replace_chain(solv_line,solvent_chain_id))
            if a not in solvent_altlocs:
                assert (a != ' '), solv_line
                solvent_altlocs.append(a) 
    solvent_lines = new_solvent_lines

    assert len(set(protein_altlocs) & set(solvent_altlocs))==len(protein_altlocs)==len(solvent_altlocs), (protein_altlocs,solvent_altlocs)
    new_altlocs={}
    for a in new_altloc_options:
        if a not in protein_altlocs:
            new_altlocs[protein_altlocs[len(new_altlocs)]]=a
        if len(new_altlocs) >= len(protein_altlocs):
            break
    else:
        raise Exception(f"{new_altlocs},{protein_altlocs}")
    def split_altloc(line):
        return new_altlocs[line[16]]


    if not sep_chain_format: # format for untangler stuff
        protein_chain_id = "A"
        for res_atom_dict in atom_dict.values():
            for d in range(2): # Split loop
                for altloc_atom_dict in res_atom_dict.values():
                    for line in altloc_atom_dict.values():
                        max_serial_num+=1
                        modified_line = replace_serial_num(line,max_serial_num)
                        modified_line = replace_chain(modified_line,protein_chain_id)
                        if d == 1:
                            modified_line = replace_altloc(modified_line,split_altloc(modified_line))
                        start_lines.append(modified_line)
    else: # Note that lines for each chain need to be contiguous in the file
        chain_dict={}
        for res_atom_dict in atom_dict.values():
            for d in range(2): # Split loop
                for altloc, altloc_atom_dict in res_atom_dict.items():
                    protein_chain_id = altloc
                    for line in altloc_atom_dict.values():
                        max_serial_num+=1
                        modified_line=line
                        modified_line = replace_serial_num(modified_line,max_serial_num)
                        modified_line = replace_chain(modified_line,protein_chain_id)
                        if d == 1:
                            altloc = split_altloc(line)
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
    for line in solvent_lines:
        solvent_resnum=int(line[22:26])
        modified_line = line
        modified_line = replace_res_num(modified_line,shift+solvent_resnum)
        modified_line=replace_serial_num(modified_line,max_serial_num)
        max_serial_num+=1
        start_lines.append(modified_line)
        splt = replace_altloc(modified_line,split_altloc(modified_line))
        splt = replace_serial_num(splt,max_serial_num)
        max_serial_num+=1
        start_lines.append(splt)

    with open(out_path,'w') as O:
        O.writelines(start_lines+end_lines)

if __name__ == "__main__":
    pdb_path = sys.argv[1]

    assert pdb_path[-4:] == ".pdb"
    out_path = pdb_path[:-4]+"_split.pdb"

    split(pdb_path,out_path)

