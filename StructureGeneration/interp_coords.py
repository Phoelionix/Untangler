#Interpolate coordinates of two models

import sys, os
from Bio.PDB import PDBParser#, parse_pdb_header
from Bio.PDB.MMCIFParser import MMCIFParser#, parse_pdb_header
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from LinearOptimizer.Input import OrderedAtomLookup


def interp_coords(out_path,modelA,modelB,weightA=0.5):
    assert 0 < weightA < 1 
    handle = os.path.basename(out_path)[:-4]
    pdbPaths = (modelA,modelB)



    builder = StructureBuilder()




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
                # switch to looking at the data in the residue of conformation i
                struct_A_residue = list(list(structs[0].get_chains())[c].get_residues())[r]
                struct_B_residue = list(list(structs[1].get_chains())[c].get_residues())[r]
                atoms={}
                for DatomA,DatomB in zip(struct_A_residue,struct_B_residue):
                    assert DatomA.get_name()==DatomB.get_name()
                    assert OrderedAtomLookup.atom_res_name(DatomA)==OrderedAtomLookup.atom_res_name(DatomB)
                    assert OrderedAtomLookup.atom_res_seq_num(DatomA)==OrderedAtomLookup.atom_res_seq_num(DatomB)
                    for atomA,atomB in zip(DatomA.__iter__(),DatomB.__iter__()):
                        assert atomA.get_altloc()==atomB.get_altloc()
                        coord = weightA*atomA.get_coord()+ (1-weightA)*atomB.get_coord()
                        bfactor = weightA*atomA.get_bfactor()+ (1-weightA)*atomB.get_bfactor()
                        occupancy = weightA*atomA.occupancy + (1-weightA)*atomB.occupancy
                        builder.init_atom(
                            atomA.get_name(),coord,bfactor,
                            occupancy=occupancy,altloc=atomA.get_altloc(),fullname=atomA.get_fullname(),
                            element=atomA.element
                            )
   

    io = PDBIO()

    builder.set_header(parse_pdb_header(pdbPaths[0]))

    io.set_structure(builder.get_structure())

    tmp_path = f"{handle}interp.pdb"
    io.save(tmp_path)

    header = ""
    with open (pdbPaths[0],'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                break
            header += line
    #print(header)

    with open(tmp_path,'r') as f:
        with open(out_path,'w') as f2: 
            f2.write(header)
            f2.write(f.read())
    os.remove(tmp_path)