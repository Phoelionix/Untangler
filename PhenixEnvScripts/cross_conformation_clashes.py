
# from cctbx.geometry_restraints import nonbonded_distance_table
# from cctbx.geometry_restraints import nonbonded_deltas
# from cctbx import geometry_restraints

# p = geometry_restraints.nonbonded_params()
# d = p.distance_table

# d.setdefault("Si")["O"] = 1.5
# print(d["Si"]["O"])


# See phenix-2.0-5793/lib/python3.9/site-packages/cctbx/geometry_restraints/tst_process_nonbonded_proxies.py

import sys,os
from cctbx.geometry_restraints import process_nonbonded_proxies,pair_proxies, manager
import mmtbx.model
import iotbx
from libtbx.utils import null_out
import cctbx.geometry_restraints.manager
from cctbx import geometry_restraints
from cctbx import uctbx, crystal
import numpy as np
import mmtbx.validation.molprobity
#from cctbx.geometry_restraints import pair_proxies
# import boost_adaptbx.boost.python as bp
# ext = bp.import_ext("cctbx_geometry_restraints_ext")
# from cctbx_geometry_restraints_ext import *

IGNORE_WATER_WATER_CLASHES=False
USE_HOLTON_CSDA=False #  csda used in scoring for untangle challenge


holton_csda = 0.6

def get_cross_conf_nonbonds(pdb_file_path,out_file,verbose,use_cdl):
    if verbose is None:
        verbose=False
    if use_cdl is None:
        use_cdl = False


    radius_considered = 10
    ###TEMPORARY####
    #radius_considered = 5
    ################
    
    params = mmtbx.model.manager.get_default_pdb_interpretation_params()
    params.pdb_interpretation.allow_polymer_cross_special_position=True
    params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = radius_considered
    params.pdb_interpretation.nonbonded_distance_cutoff= radius_considered
    params.pdb_interpretation.restraints_library.cdl = use_cdl
    if USE_HOLTON_CSDA:
        params.pdb_interpretation.const_shrink_donor_acceptor=holton_csda
    #pdb_inp = iotbx.pdb.input(lines=raw_records.split("\n"), source_info=None)

    tmp_pdb_file = os.path.join(os.path.abspath(os.path.join(__file__ ,"../")),"tmp_samealtloc.pdb")

    first_altloc={}
    with open(pdb_file_path) as f, open(tmp_pdb_file,"w") as w:
        lines = f.readlines()
        max_resnum=num_altlocs=0
        conformation_number:dict[str,int]={}
        for line in lines:
            valid_record_types=["ATOM","HETATM"]
            if any([line.startswith(k) for k in valid_record_types]):
                og_resnum=int(line[22:26])
                max_resnum= max(max_resnum,og_resnum)
                altloc = line[16]
                
                if og_resnum not in first_altloc:
                    first_altloc[og_resnum]=altloc
        for line in lines:
            valid_record_types=["ATOM","HETATM"]
            if any([line.startswith(k) for k in valid_record_types]):
                og_resnum = int(line[22:26])
                altloc = line[16]
                #new_chain="X"
                new_chain=old_chain=line[21]
                if altloc !=first_altloc[og_resnum]:
                    continue
                new_resnum=og_resnum

                if line[17:20]=="HOH":
                    line= "HETATM"+line[6:]
                w.write(line[:16]+"A"+line[17:21]+new_chain+ f"{new_resnum}".rjust(4)+line[26:])
            else:
                w.write(line)
            #same_altloc_labels.append(flex.std_string([new_str]))
    single_conformation_pdb_inp = iotbx.pdb.input(tmp_pdb_file)


    #pdb_inp = iotbx.pdb.input(pdb_file_path)


    model = mmtbx.model.manager(
        model_input = single_conformation_pdb_inp,
        log         = null_out(),
        # crystal_symmetry=crystal.symmetry(
        #         unit_cell=unreasonably_large_unit_cell,
        #         space_group_symbol="P1")
    )
    model.process(pdb_interpretation_params=params,
        make_restraints=True)


    


    #pnp_manager = process_nonbonded_proxies.manager(model=model)
    # grm = pnp_manager.model.get_restraints_manager().geometry
    # xrs = pnp_manager.model.get_xray_structure()
    # sites_cart  = pnp_manager.model.get_sites_cart()

    grm: cctbx.geometry_restraints.manager.manager
    grm = model.get_restraints_manager().geometry
    #grm = validation.model.get_restraints_manager().geometry
    xrs = model.get_xray_structure()
    sites_cart  = model.get_sites_cart()
    site_labels = xrs.scatterers().extract_labels()





    pair_proxies: geometry_restraints.pair_proxies
    

    #grm.crystal_symmetry=None
    pair_proxies = grm.pair_proxies(
                        sites_cart  = sites_cart,
                        site_labels = site_labels)
    proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted( # returns C++ nonbonded_sorted_asu_proxies
        by_value    = "delta",
        sites_cart  = sites_cart,
        site_labels = site_labels)
    #print(proxies_info_nonbonded[0][0])


    assert proxies_info_nonbonded is not None
    nonbonded_list = proxies_info_nonbonded[0]



    og_site_labels= mmtbx.model.manager(
        model_input = iotbx.pdb.input(pdb_file_path),
        log         = null_out()
    ).get_xray_structure().scatterers().extract_labels()
    og_sites_cart  = mmtbx.model.manager(
        model_input = iotbx.pdb.input(pdb_file_path),
        log         = null_out()
    ).get_sites_cart()
    #og_site_labels=site_labels


    # consider all altlocs (lazy implementation)
    ordered_atom_sites_dict={}
    def sep(X,Y):
        return np.sqrt(np.sum((np.array(X)-np.array(Y))**2))
    for site_label,xyz in zip(og_site_labels,og_sites_cart):
        key = site_label[:9]+site_label[10:]
        
        if key not in ordered_atom_sites_dict:
            ordered_atom_sites_dict[key]=[]
        ordered_atom_sites_dict[key].append((site_label,xyz))
        #print(key,site_label,xyz)

        
    def read_clash_table(keyA,keyB): #  https://doi.org/10.1002/pro.3330
        vdw_dict=dict(
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

    out_data=[]
    #print("WARNING clash table is not  fully implemented")
    for i, item in enumerate(nonbonded_list):
        if i%100000==0:
            print(f"{i}/{len(nonbonded_list)}")
        i_seq          = item[1]
        j_seq          = item[2]
        #model_distance = item[3]
        #vdw_sum        = item[4]
        symop_str      = item[5] 
        #symop          = item[6]




        keyA,keyB = [site_label[:9]+site_label[10:] for site_label in [site_labels[i_seq],site_labels[j_seq]]]
        # if keyA=='pdb=" C  VAL A   1 "':
        #     if keyB[-4:-2]==" 2":
        #         print(keyA,keyB)
        #         print(site_labels[i_seq],site_labels[j_seq])
        
        vdw_sum=read_clash_table(keyA,keyB)
        if vdw_sum is None:
            continue

        debug_packing_added={}
        debug_nonpacking_added={}
        conformer_site_label_A,conformer_site_label_B = ordered_atom_sites_dict[keyA][0][0],ordered_atom_sites_dict[keyB][0][0]

        is_crystal_packing_contact= "true" if symop_str.strip()!="" else "false"

        datum=[
            conformer_site_label_A.split('"')[1], #pdb1 (pdb entry 1)
            conformer_site_label_B.split('"')[1], #pdb2 (pdb entry 2)
            vdw_sum,
            is_crystal_packing_contact,
        ]
        debug_key=datum[0]+" " +datum[1]
        debug_dict = debug_packing_added if is_crystal_packing_contact else debug_nonpacking_added
        # Assuming only difference in vdw for same energy types is whether it is a crystal packing contact
        if debug_key in debug_dict:
            assert debug_dict[debug_key]==vdw_sum, (datum,debug_dict[debug_key])
        else:
            debug_dict[debug_key]=vdw_sum
            out_data.append(datum)
            
    if out_file is not None:
        with open(out_file,"w") as f:
            f.write('\n'.join(['|'.join([str(i) for i in items]) for items in out_data]))
    print(f"number of cross-conformer clash entries: {len(out_data)}")
    if verbose:
        pair_proxies.nonbonded_proxies.show_histogram_of_model_distances(
            sites_cart=sites_cart,
            f=sys.stdout)
    return out_data
if __name__=="__main__":
    assert len(sys.argv)>=2
    args = [sys.argv[1],None,None,None]
    for i, arg in enumerate(sys.argv[1:]):
        args[i]=arg
    get_cross_conf_nonbonds(*args)