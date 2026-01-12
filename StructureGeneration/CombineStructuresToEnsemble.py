
import sys, os
from Bio.PDB import PDBParser#, parse_pdb_header
from Bio.PDB.MMCIFParser import MMCIFParser#, parse_pdb_header
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder

handle = sys.argv[1] 
pdbPaths = sys.argv[2:]



builder = StructureBuilder()


altlocs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(pdbPaths)]


builder.init_structure("ensemble")
builder.init_model("M")



# first_pdb_path = pdbPaths[0]
# other_pdb_paths = pdbPaths[1:]

structs = [PDBParser().get_structure("struct",pdb_path)  for pdb_path in pdbPaths]

# structs = []
# for pdb_path in pdbPaths:
#     if pdb_path[-3:]==".pdb":
#         parser=PDBParser()
#     if pdb_path[-3:]==".cif":
#         parser=MMCIFParser()
#     else:
#         assert False
#     structs.append(parser.get_structure("struct",pdb_path))

ignoreHetero=False
for c, chain in enumerate(structs[0].get_chains()):  # loop through just the first struct for chain and residue, we only want to have differences for atoms
    builder.init_chain(chain.get_id())
    for r, res in enumerate(chain.get_residues()): 
        builder.init_seg(res.get_segid())
        builder.init_residue(res.get_resname(),*res.get_id())
        if res.get_id()[0] == " " or not ignoreHetero: 
            for i, altloc in enumerate(altlocs):
                # switch to looking at the data in the residue of conformation i
                conformation_residue = list(list(structs[i].get_chains())[c].get_residues())[r]
                for atom in conformation_residue.get_atoms():
                    builder.init_atom(
                    atom.get_name(),atom.get_coord(),atom.get_bfactor(),
                    occupancy=atom.occupancy/len(altlocs),altloc=altloc,fullname=atom.get_fullname(),
                    element=atom.element
                    )

io = PDBIO()

builder.set_header(parse_pdb_header(pdbPaths[0]))

io.set_structure(builder.get_structure())

out_path = f"{handle}.pdb"
tmp_path = "structure_combined_tmp.pdb"
io.save(tmp_path)



header = ""
with open (pdbPaths[0],'r') as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            break
        header += line
print(header)

with open(tmp_path,'r') as f:
    with open(out_path,'w') as f2: 
        f2.write(header)
        f2.write(f.read())
os.remove(tmp_path)