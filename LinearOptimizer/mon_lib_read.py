# https://www.ccp4.ac.uk/html/mon_lib.html

import subprocess
import os
from UntangleFunctions import UNTANGLER_WORKING_DIRECTORY


#TODO switch to vdw radius?
def read_lj_parameters():

    params:dict[str,dict[str,tuple[float]]]={}
    with open(f"{UNTANGLER_WORKING_DIRECTORY}../ccp4-8.0/lib/data/monomers/ener_lib.cif") as f:
        print("lj param warning: Reading parameters for atoms with hydrogens")
        reading=False
        for line in f:
            if line.startswith("#             ENERGY =  EPSij * ( (Rmin/Rij)**12 - 2 * (Rmin/Rij)**6 )"):
                reading=True
            elif line.startswith("# ---  H-BONDS ---"):
                reading=False
                break
            elif line.startswith("#"):
                continue
            elif reading:
                A,B,E_min,R_min,h_flag = line.split()
                E_min,R_min = float(E_min),float(R_min)
                for a, b in ((A,B),(B,A)):
                    if a not in params:
                        params[a]={}
                    # H flag params have priority. (Assume has H where possible for now).
                    if b in params[a] and h_flag!= "h":
                        continue  
                    params[a][b]=(E_min,R_min)
        assert not reading
    return params
    
def read_vdw_radii(with_H):

    params:dict[str,dict[float]]={}
    with open(f"{UNTANGLER_WORKING_DIRECTORY}../ccp4-8.0/lib/data/monomers/ener_lib.cif") as f:
        reading=False
        for line in f:
            if line.startswith("#   _vdw_radius  Van-Der-Waals radius"):
                reading=True
            elif line.startswith("#   Sort out valencies of atoms"):
                reading=False
                break
            elif line.startswith("#"):
                continue
            elif reading:
                A,_,_,r_vdw,r_vdw_H,_,_,_,_ = line.split()
                if with_H and r_vdw_H.isnumeric():
                    params[A]=float(r_vdw_H)
                else:
                    params[A]=float(r_vdw)
        assert not reading
    return params

def get_mon_lib_names(residue):
    # Returns a dict that gives the atom type used in ener_lib.cif for each corresponding pdb atom name (as defined by MON_LIB)
    out_folder = os.path.join(UNTANGLER_WORKING_DIRECTORY,"LinearOptimizer","tmp_out","")
    libcheck_script = os.path.join(UNTANGLER_WORKING_DIRECTORY,"LinearOptimizer","mon_libcheck.sh")
    residue_monomer_lib_file= os.path.join(out_folder,residue+".lib")
    if not os.path.exists(residue_monomer_lib_file):
        args = ["bash",libcheck_script,residue,out_folder] 
        print (f"|+ Running: {' '.join(args)}")
        subprocess.run(args)#,stdout=log)

    vdw_dict={}
    lj_dict={}
    with open(residue_monomer_lib_file) as f:
        reading=False
        lj_substitutes = dict(CR6="C",CR5="C",CR15="CH1",CR56="C",NT1="NH1",NR15="NH1",NR16="NH1",NC1="NH1",NC2="NH2",NC3="NH3",NH3="NT3") # Atom energy types that dont appear in vdw contacts (for some probably because charged... but doing this for now in lieu of Coulomb potential...)???  (NB NH3/NT3 is synonymous)
        for line in f:
            if line.startswith("_chem_comp_atom.z"):
                reading=True
            elif reading:
                if line.startswith("loop_"):
                    reading=False
                    break
                name = line.split()[1]
                type_energy=line.split()[3]
                vdw_dict[name]=type_energy
                if type_energy in lj_substitutes:
                    type_energy=lj_substitutes[type_energy]
                lj_dict[name]=type_energy                
        assert not reading
    return vdw_dict,lj_dict

            # 4th element: _chem_comp_atom.type_energy
            





# _lib_vdw.atom_type_1 atomic chemical type
# _lib_vdw.atom_type_2
# _lib_vdw.energy_min EPSij  minimum of energy parameter
# _lib_vdw.radius_min Rmin radius of the minimum of energy parameter
# _lib_vdw.H_flag "h" - the parameters for atoms with hydrogens 
#                  Rij - actual distance
#  _H_flag         "h" - the parameters for atoms with hydrogens 
#

#             ENERGY =  EPSij * ( (Rmin/Rij)**12 - 2 * (Rmin/Rij)**6 )
#

#C    C       -0.19686     3.600   .