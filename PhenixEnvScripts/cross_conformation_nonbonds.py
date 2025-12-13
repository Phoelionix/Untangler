
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
#from cctbx.geometry_restraints import pair_proxies
# import boost_adaptbx.boost.python as bp
# ext = bp.import_ext("cctbx_geometry_restraints_ext")
# from cctbx_geometry_restraints_ext import *


holton_csda = 0.6

def get_cross_conf_nonbonds(pdb_file_path,out_file,verbose,use_cdl):
    if verbose is None:
        verbose=False
    if use_cdl is None:
        use_cdl = False

    params = mmtbx.model.manager.get_default_pdb_interpretation_params()
    params.pdb_interpretation.allow_polymer_cross_special_position=True
    params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = 15
    params.pdb_interpretation.nonbonded_distance_cutoff= 15
    params.pdb_interpretation.restraints_library.cdl = use_cdl
    params.pdb_interpretation.const_shrink_donor_acceptor=holton_csda
    #pdb_inp = iotbx.pdb.input(lines=raw_records.split("\n"), source_info=None)

    tmp_pdb_file = os.path.join(os.path.abspath(os.path.join(__file__ ,"../")),"tmp_samealtloc.pdb")

    with open(pdb_file_path) as f, open(tmp_pdb_file,"w") as w:
        lines = f.readlines()
        max_resnum=num_altlocs=0
        conformation_number:dict[str,int]={}
        for line in lines:
            valid_record_types=["ATOM","HETATM"]
            if any([line.startswith(k) for k in valid_record_types]):
                max_resnum= max(max_resnum,int(line[22:26]))
                altloc = line[16]
                if len(conformation_number)>0:  # single altloc
                    continue
                if altloc not in conformation_number:
                    num_altlocs+=1
                    conformation_number[altloc]=num_altlocs
        for line in lines:
            valid_record_types=["ATOM","HETATM"]
            if any([line.startswith(k) for k in valid_record_types]):
                resnum = int(line[22:26])
                altloc = line[16]
                #new_chain="X"
                new_chain=old_chain=line[21]
                if altloc not in conformation_number:
                    continue
                new_resnum=resnum+(-1+conformation_number[altloc])*max_resnum
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

        

    out_data=[]
    for i, item in enumerate(nonbonded_list):
        if i%100000==0:
            print(f"{i}/{len(nonbonded_list)}")
        i_seq          = item[1]
        j_seq          = item[2]
        #model_distance = item[3]
        vdw_sum        = item[4]
        symop_str      = item[5] 
        symop          = item[6]


        keyA,keyB = [site_label[:9]+site_label[10:] for site_label in [site_labels[i_seq],site_labels[j_seq]]]
        # if keyA=='pdb=" C  VAL A   1 "':
        #     if keyB[-4:-2]==" 2":
        #         print(keyA,keyB)
        #         print(site_labels[i_seq],site_labels[j_seq])
        distance_cutoff = 4
        
        debug_packing_added={}
        debug_nonpacking_added={}
        for (conformer_site_label_A,coordA) in ordered_atom_sites_dict[keyA]:
            for (conformer_site_label_B,coordB) in ordered_atom_sites_dict[keyB]:
                model_distance=sep(coordA,coordB)
                if sep(coordA,coordB) > distance_cutoff:
                    continue
                # if keyA=='pdb=" C  GLY A  59 "':
                #     if keyB[-5:-2]==" 60":
                #         print(model_distance,vdw_sum)
                #         print(conformer_site_label_A,conformer_site_label_B)
                #         print(site_labels[i_seq],site_labels[j_seq])
                crystal_packing_contact= "true" if symop_str.strip()!="" else "false"
    
                datum=[
                    conformer_site_label_A.split('"')[1], #pdb1 (pdb entry 1)
                    conformer_site_label_B.split('"')[1], #pdb2 (pdb entry 2)
                    vdw_sum,
                    crystal_packing_contact,
                ]
                debug_key=datum[0]+" " +datum[1]
                debug_dict = debug_packing_added if crystal_packing_contact else debug_nonpacking_added
                # Assuming only difference in vdw for same energy types is whether it is a crystal packing contact
                if debug_key in debug_dict:
                    assert debug_dict[debug_key]==vdw_sum, (datum,debug_dict[debug_key])
                else:
                    debug_dict[debug_key]=vdw_sum
                    out_data.append(datum)
            
    if out_file is not None:
        with open(out_file,"w") as f:
            f.write('\n'.join(['|'.join([str(i) for i in items]) for items in out_data]))
    print(f"number of cross-conformer nonbonds: {len(out_data)}")
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