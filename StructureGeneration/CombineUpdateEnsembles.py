
import sys, os
from Bio.PDB import PDBParser#, parse_pdb_header
from Bio.PDB.MMCIFParser import MMCIFParser#, parse_pdb_header
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Atom import Atom,DisorderedAtom

out_path = sys.argv[1] 
pdbPaths = sys.argv[2::2]  
altlocs_to_use = sys.argv[3::2]

if len(sys.argv) < 6:
    print("Usage: python3.9 path/to/output.pdb path/to/pdb1.pdb AC path/to/pdb2.pdb BDEFGH")



# print(pdbPaths)
# print(altlocs_to_use)

for i, a in enumerate(altlocs_to_use):
    for j, b in enumerate(altlocs_to_use):
        if j <=i: continue
        assert len(set(a)&set(b))==0

# Only outputs based on self.better_dict at present!!
reference_pdb_file=pdbPaths[0]
with open (out_path,'w') as O:
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
            O.write(line)
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
                        done=True
                        break
                else: # Not an atom entry
                    continue
                # is an atom entry
                altloc = line[16]
                if altloc in allowed_altlocs:
                    O.write(replace_serial_num(line,next_sn))
                    next_sn+=1

    O.write("END\n")
