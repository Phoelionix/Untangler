
import sys, os
from LinearOptimizer.Solver import MTSP_Solver
import untangle,UntangleFunctions



def get_altlocs(pdb_path):         
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
            if resname in solvent_res_names:
                if altloc not in solvent_altlocs:
                    solvent_altlocs.append(altloc) 
                continue
                
            if altloc not in protein_altlocs:
                protein_altlocs.append(altloc) 
    return protein_altlocs,solvent_altlocs


def split_worst(pdb_path,out_path,split_waters,sep_chain_format=False,protein_altloc_from_chain_fix=False,missing_water_altloc_fix=True):
    assert False, "Bad idea, want interchangeable occupancies"

    altloc_options = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    protein_altlocs, water_altlocs = get_altlocs(pdb_path)

    new_altloc=None
    for a in (altloc_options):
        if a not in (protein_altlocs+water_altlocs):
            new_altloc = a
            break
    else:
        assert False, "No available new altloc found"

    MTSP_Solver.prepare_geom_files(pdb_path,protein_altlocs,water_swaps=False)
    worst_score = 0
    worst_altloc = None
    for altloc in protein_altlocs:
        subset_model = MTSP_Solver.subset_model_path(pdb_path,altloc)
        score = untangle.Untangler.Score(*UntangleFunctions.get_score(UntangleFunctions.score_file_name(subset_model)))
        print(f"Altloc: {altloc}, Score: {score}")
        if score.combined >= worst_score:
            worst_score = score.combined
            worst_altloc = altloc

    print(f"Worst scoring altloc: {worst_altloc} ({worst_score})")

    
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

    new_altlocs={worst_altloc: new_altloc}

    def split_altloc(line):
        original_altloc = line[16]
        return new_altlocs[original_altloc] if original_altloc in new_altlocs else None 


    if not sep_chain_format: # format for untangler stuff
        protein_chain_id = "A"
        for res_atom_dict in atom_dict.values():
            for d in range(2): # Split loop
                for altloc_atom_dict in res_atom_dict.values():
                    for line in altloc_atom_dict.values():
                        max_serial_num+=1
                        modified_line = replace_serial_num(line,max_serial_num)
                        modified_line = replace_chain(modified_line,protein_chain_id)
                        child_altloc = split_altloc(modified_line)
                        if child_altloc is not None:
                            # Split
                            modified_line = replace_occupancy(line,occupancy/2)

                        if d == 1:
                            child_altloc = split_altloc(modified_line)
                            if child_altloc is None: #  i.e. if we aren't splitting this entry into two altlocs
                                continue
                            modified_line = replace_altloc(modified_line,child_altloc)
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
                            child_altloc = split_altloc(modified_line)
                            if child_altloc is None:  #  i.e. if we aren't splitting this entry into two altlocs
                                continue
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
        
        child_altloc = split_altloc(modified_line)
        do_split = (child_altloc is not None and split_waters)
        if do_split:
            occupancy=float(modified_line[54:60])
            modified_line = replace_occupancy(modified_line,occupancy/2)

        modified_line = replace_res_num(modified_line,shift+solvent_resnum)
        modified_line=replace_serial_num(modified_line,max_serial_num)
        max_serial_num+=1
        start_lines.append(modified_line)

        if do_split:
            splt = replace_altloc(modified_line,child_altloc)
            splt = replace_serial_num(splt,max_serial_num)
            max_serial_num+=1
            start_lines.append(splt)

    with open(out_path,'w') as O:
        O.writelines(start_lines+end_lines)

if __name__ == "__main__":
    pdb_path = sys.argv[1]

    assert pdb_path[-4:] == ".pdb"
    out_path = pdb_path[:-4]+"_splitworst.pdb"

    split_worst(pdb_path,out_path,split_waters=False)

