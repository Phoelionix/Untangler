IGNORE_WATER_WATER_CLASHES=False

VDW_RADII=dict( #  https://doi.org/10.1002/pro.3330
    H=1.22, H_aromatic=1.05, H_polar=1.05,
    C=1.7,
    C_aromatic=1.75,
    N=1.55,
    O=1.40, HOH_O=1.40,
    P=1.80,
    S=1.80,
    MN_MN=1.80, # FIXME
    #Se=1.90,
)

def pdb_key(atom_name:str,resname:str,chain:str="A",resnum:int=1)->str:
    return f'pdb="{atom_name:<4}{resname.strip():>3} {chain}{resnum:>4} "'

def read_clash_table(keyA,keyB): 
    vdw_dict=VDW_RADII
    #TODO HOH
    # TODO aromatic, polar
    # key of form: 'pdb=" C  VAL A   1 "'


    def to_energy_type(key):
        entry=key.split('"')[1]
        name,resname = entry[0:4].strip(),entry[4:7].strip()

        aromatic_H_dict=dict(
            HIS=("HD1","HE1","HE2","HD2"),
            TYR=("HD1","HE1","HE2","HD2"),
            PHE=("HD1","HE1","HE2","HD2","HZ"),
            TRP=("HD1","HE1","HH2","HZ2","HZ3","HE3"),
            NAP=("H8A","H2A","H2N","H4N","H5N","H6N"),
            FOL=("HN1","H7","H15","H12","H13","H15","H16"),
        )
        aromatic_C_dict=dict(
            HIS=("CG","CD2","CE1"),
            TYR=("CG","CD1","CD2","CE1","CE2","CZ"),
            PHE=("CG","CD1","CD2","CE1","CE2","CZ"),
            TRP=("CD1","CG","CD2","CE2","CZ2","CH2","CZ3","CE3"),
            NAP=("C2A","C4A","C5A","C6A","C8A","C2N","C3N","C4N","C5N","C6N"),
            FOL=("C2","C4","C4A","C8A","C6","C7","C11","C12","C13","C14","C15","C16"),
        )
        polar_H_dict=dict(
            ARG=("HE","HH11","HH12","HH21","HH22"),
            LYS=("HZ1","HZ2","HZ3"),
            SER=("HG",),
            GLN=("HE21","HE22"),
            TRP=("HE1"),
            CSD=("HD2"),
            HIS=("HD1","HE2"),
            TYR=("HH"),
            ASN=("HD21","HD22"),
            NAP=("H71N","H72N","HO2N","HO3A","H61A","H62A"),
            FOL=("HN21","HN22","HN1","HN0","HN","HOE2",),
            THR=("HG1",),
            CYS=("HG",),
        )
        if resname in polar_H_dict and  name in polar_H_dict[resname]\
            or name in ["H","HN","H1","H2","H3"]:
            return "H_polar"
        if resname in aromatic_H_dict and name in aromatic_H_dict[resname]:
            return "H_aromatic"
        if resname in aromatic_C_dict and name in aromatic_C_dict[resname]:
            return "C_aromatic"

        return f"{resname}_{name}" if f"{resname}_{name}" in vdw_dict else name[0]
    def is_acceptor(key):
        # Oxygen, or nitrogen not bonded to H
        entry=key.split('"')[1]
        name,resname = entry[0:4].strip(),entry[4:7]
        if name[0]=="O":
            return True
        if resname == "NAP" and name in ("N3A","N9A","N7A","N1A","N1N"):
            return True
        if resname == "FOL" and name in ("N3","N5","N8"):
            return True



    resnameA,resnameB = keyA.split('"')[1][4:7],keyB.split('"')[1][4:7]
    if IGNORE_WATER_WATER_CLASHES:
        if resnameA==resnameB=="HOH":
            return None


    def possible_hydrogen_bond(keyA,keyB):
        return "H_polar" in (to_energy_type(keyA),to_energy_type(keyB))\
            and  (is_acceptor(keyA) or is_acceptor(keyB))

    if possible_hydrogen_bond(keyA,keyB):    
        return None
    try: 
        return vdw_dict[to_energy_type(keyA)]+vdw_dict[to_energy_type(keyB)]
    except:
        return None
