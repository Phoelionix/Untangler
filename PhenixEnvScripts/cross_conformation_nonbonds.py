
# from cctbx.geometry_restraints import nonbonded_distance_table
# from cctbx.geometry_restraints import nonbonded_deltas
# from cctbx import geometry_restraints

# p = geometry_restraints.nonbonded_params()
# d = p.distance_table

# d.setdefault("Si")["O"] = 1.5
# print(d["Si"]["O"])


# See phenix-2.0-5793/lib/python3.9/site-packages/cctbx/geometry_restraints/tst_process_nonbonded_proxies.py

# if not running phenix.python... and just python3.9...

import sys,os
# # old_sys_path = [p for p in sys.path]

# # sys.path=[os.path.join(os.environ["PHENIX"],"lib","python3.9","site-packages")]


# phenix_lib = os.path.join(os.environ["PHENIX"],"lib","python3.9","site-packages")

# if phenix_lib not in sys.path:
#     sys.path.append(phenix_lib)

# # Ensure required Phenix environment variables exist
# #os.environ.setdefault("PHENIX", phenix_root)
# os.environ.setdefault("LIBTBX_BUILD", os.path.join(os.environ["PHENIX"], "share","cctbx"))
# #os.environ.setdefault("CCTBX_BUILD", os.path.join(phenix_root, "share"))
# #os.environ["PHENIX_RESOURCES"] = os.path.join(phenix_root, "share")

from cctbx.geometry_restraints import process_nonbonded_proxies,pair_proxies
import mmtbx.model
import iotbx
from libtbx.utils import null_out




def get_cross_conf_nonbonds(pdb_file_path,out_file=None,verbose=False):


    params = mmtbx.model.manager.get_default_pdb_interpretation_params()
    params.pdb_interpretation.allow_polymer_cross_special_position=True
    params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
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
                if altloc not in conformation_number:
                    num_altlocs+=1
                    conformation_number[altloc]=num_altlocs
        for line in lines:
            valid_record_types=["ATOM","HETATM"]
            if any([line.startswith(k) for k in valid_record_types]):
                resnum = int(line[22:26])
                altloc = line[16]
                new_resnum=resnum+(-1+conformation_number[altloc])*max_resnum
                w.write(line[:16]+"x"+line[17:22]+ f"{new_resnum}".rjust(4)+line[26:])
            else:
                w.write(line)
            #same_altloc_labels.append(flex.std_string([new_str]))


    pdb_inp = iotbx.pdb.input(tmp_pdb_file)
    model = mmtbx.model.manager(
        model_input = pdb_inp,
        log         = null_out()
    )
    model.process(pdb_interpretation_params=params,
        make_restraints=True)





    #pnp_manager = process_nonbonded_proxies.manager(model=model)
    # grm = pnp_manager.model.get_restraints_manager().geometry
    # xrs = pnp_manager.model.get_xray_structure()
    # sites_cart  = pnp_manager.model.get_sites_cart()

    grm = model.get_restraints_manager().geometry
    xrs = model.get_xray_structure()
    sites_cart  = model.get_sites_cart()
    site_labels = xrs.scatterers().extract_labels()






    pair_proxies = grm.pair_proxies(
                        sites_cart  = sites_cart,
                        site_labels = site_labels)
    proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted(
        by_value    = "delta",
        sites_cart  = sites_cart,
        site_labels = site_labels)


    assert proxies_info_nonbonded is not None
    nonbonded_list = proxies_info_nonbonded[0]



    og_site_labels= mmtbx.model.manager(
        model_input = iotbx.pdb.input(pdb_file_path),
        log         = null_out()
    ).get_xray_structure().scatterers().extract_labels()

    num_skipped=0
    out_data=[]
    for item in nonbonded_list:
        i_seq          = item[1]
        j_seq          = item[2]
        model_distance = item[3]
        vdw_sum        = item[4]
        symop_str      = item[5] 
        symop          = item[6]

        # Skip if same atom but different altloc
        if og_site_labels[i_seq][:9]+og_site_labels[i_seq][10:]==og_site_labels[j_seq][:9]+og_site_labels[j_seq][10:]:
            num_skipped+=1    
            continue
        out_data.append([
            og_site_labels[i_seq].split('"')[1], #pdb1 (pdb entry 1)
            og_site_labels[j_seq].split('"')[1], #pdb2 (pdb entry 2)
            vdw_sum,
        ])
        if verbose:
            print(og_site_labels[i_seq],og_site_labels[j_seq],float(vdw_sum))
    if out_file is not None:
        with open(out_file,"w") as f:
            f.write('\n'.join(['|'.join([str(i) for i in items]) for items in out_data]))
    print(f"number of cross-conformer nonbonds: {len(nonbonded_list)-num_skipped}")
    if verbose:
        pair_proxies.nonbonded_proxies.show_histogram_of_model_distances(
            sites_cart=sites_cart,
            f=sys.stdout)
    return out_data
if __name__=="__main__":
    assert len(sys.argv)>=2
    args = [sys.argv[1],None,None]
    for i, arg in enumerate(sys.argv[1:]):
        args[i]=arg
    get_cross_conf_nonbonds(*args)

# That was nearly instant.
# Now need to just modify phenix to generate vdw sums for all possible altlocs or find the function that does

