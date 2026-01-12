
import sys, os
from Bio.PDB import PDBParser#, parse_pdb_header
from Bio.PDB.MMCIFParser import MMCIFParser#, parse_pdb_header
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Atom import Atom,DisorderedAtom





def combine_update_ensembles(out_path,pdbPaths,altlocs_to_use:list[str]):
    # if ann element of altlocs_to_use is None, it means use all altlocs that aren't specified anywhere in the list.
    
    # print("=============")
    # print(pdbPaths)
    # print(altlocs_to_use)
    # print(out_path)
    # print("=============")
    assert len(pdbPaths)==len(altlocs_to_use)
    for i, a in enumerate(altlocs_to_use):
        if a is None: continue
        for j, b in enumerate(altlocs_to_use):
            if b is None: continue
            if j <=i: continue
            assert len(set(a)&set(b))==0

    wild_idx=None
    nonwild_altlocs=""
    
    
    for i, allowed_altlocs in enumerate(altlocs_to_use):            
        if allowed_altlocs == [None] or allowed_altlocs is None:
            altlocs_to_use[i]=[None]
            if wild_idx is not None:
                raise Exception("No more than 1 altloc selection can be left wild")
            wild_idx=i
        else:
            nonwild_altlocs+=allowed_altlocs


    # Only outputs based on self.better_dict at present!!
    out_str=""
    reference_pdb_file=pdbPaths[0]
    atom_line_start_strings = ["ATOM","HETATM"]
    with open(reference_pdb_file) as R:
        for line in R:
            done=False
            for s in atom_line_start_strings:
                if line.startswith(s):
                    done=True
                    break
            if done:
                break
            out_str+=line
            continue
    
    def replace_serial_num(line,serial_num):
        serial_num = str(serial_num)
        serial_num = ' '*(5-len(serial_num))+serial_num
        return line[:6]+serial_num+line[11:]
    next_sn=1
    for pdb_file,allowed_altlocs in zip(pdbPaths,altlocs_to_use):
        with open(pdb_file) as R:        
            for line in R:
                for s in atom_line_start_strings:
                    if line.startswith(s):
                        break
                else: 
                    # Not an atom entry
                    continue
                # is an atom entry
                altloc = line[16]
                assert altloc is not None
                wrote=False
                if (altloc in allowed_altlocs) or (allowed_altlocs == [None] and altloc not in nonwild_altlocs):
                    out_str+= replace_serial_num(line,next_sn)
                    next_sn+=1
                    wrote=True

    out_str+="END\n"
    with open (out_path,'w') as O:
        O.write(out_str)

if __name__ == "__main__":

    out_path = sys.argv[1] 
    pdbPaths = sys.argv[2::2]  
    altlocs_to_use = sys.argv[3::2]

    if len(sys.argv) < 6:
        print("Usage: python3.9 path/to/output.pdb path/to/pdb1.pdb AC path/to/pdb2.pdb BDEFGH")
    combine_update_ensembles(out_path,pdbPaths,altlocs_to_use)