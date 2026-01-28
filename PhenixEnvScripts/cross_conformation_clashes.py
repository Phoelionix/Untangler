
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
from mmtbx.validation.clashscore import clashscore
from libtbx import easy_run
import libtbx
import json
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

    SINGLE_ALTLOC_MODE=False
    new_resnum_to_old_resnum_dict={}
    with open(pdb_file_path) as f, open(tmp_pdb_file,"w") as w:
        lines = f.readlines()
        max_resnum=num_altlocs=0
        conformation_number:dict[str,int]={}
        for line in lines:
            valid_record_types=["ATOM","HETATM"]
            if any([line.startswith(k) for k in valid_record_types]):
                max_resnum= max(max_resnum,int(line[22:26]))
                altloc = line[16]
                # if altloc != "B":
                #     continue
                if len(conformation_number)>0 and SINGLE_ALTLOC_MODE:  # single altloc
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
                #new_resnum=resnum+(-1+conformation_number[altloc])*max_resnum
                new_resnum=resnum+(-1+conformation_number[altloc])*max_resnum
                new_resnum_to_old_resnum_dict[new_resnum]=resnum

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
    print("Processing")
    model.process(pdb_interpretation_params=params,
        make_restraints=True)
    print("Processed")

    


    #pnp_manager = process_nonbonded_proxies.manager(model=model)
    # grm = pnp_manager.model.get_restraints_manager().geometry
    # xrs = pnp_manager.model.get_xray_structure()
    # sites_cart  = pnp_manager.model.get_sites_cart()

    #grm: cctbx.geometry_restraints.manager.manager
    #grm = model.get_restraints_manager().geometry
    #grm = validation.model.get_restraints_manager().geometry
    xrs = model.get_xray_structure()
    sites_cart  = model.get_sites_cart()
    site_labels = xrs.scatterers().extract_labels()





    #pair_proxies: geometry_restraints.pair_proxies
    

    #grm.crystal_symmetry=None
    #validation =  mmtbx.validation.molprobity.molprobity(model=model)
    # grm = validation.model.get_restraints_manager().geometry
    # pair_proxies = grm.pair_proxies(
    #                     sites_cart  = sites_cart,
    #                     site_labels = site_labels)
    # proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted( # returns C++ nonbonded_sorted_asu_proxies
    #     by_value    = "delta",
    #     sites_cart  = sites_cart,
    #     site_labels = site_labels)
    #print(proxies_info_nonbonded[0][0])
    #validation.show_summary()



    CLASH_SCAN_BUFFER=0 # Doesnt seem to work
    def run_probe_clashscore(probe_clashscore_manager, pdb_string):
        probe_cmd=libtbx.env.under_build(os.path.join('probe', 'exe', 'probe'))

        add_vdw_arg = f"-ADDvdw{int(CLASH_SCAN_BUFFER)}.0"
        #probe_cmd_and_args=f'{probe_cmd} -u -q -mc -het -once {add_vdw_arg}  "ogt0 not water" "ogt0" - ' 
        probe_cmd_and_args=f'{probe_cmd} -u -q -mc -het -once {add_vdw_arg}  "ogt0 not water" "ogt0" - ' 

        print("Running:",probe_cmd_and_args, "[pdb stdin_lines]")

        probe_out = easy_run.fully_buffered(probe_cmd_and_args,stdin_lines=pdb_string)
        if (probe_out.return_code != 0):
            raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
            "\n".join(probe_out.stderr_lines))
        probe_unformatted = probe_out.stdout_lines

        filtered_dicts_list=probe_clashscore_manager.process_raw_probe_output(probe_unformatted)
        return filtered_dicts_list




    geo = mmtbx.model.statistics.geometry(
      model           = model,
      fast_clash      = False,
      condensed_probe = False)
    clash_score = clashscore(
        pdb_hierarchy   = geo.pdb_hierarchy,
        condensed_probe = False,
        fast            = False,
        keep_hydrogens  = geo.use_hydrogens,
        nuclear         = geo.model.is_neutron())
    pdb_string=clash_score.probe_clashscore_manager.h_pdb_string
    #print(clash_score.probe_clashscore_manager.full_probe_txt, pdb_string)
    filtered_dicts = run_probe_clashscore(clash_score.probe_clashscore_manager,pdb_string)

    #print(dir(filtered_dicts[0]))


    atom_xyz_dict={}
    def sep(X,Y):
        return np.sqrt(np.sum((np.array(X)-np.array(Y))**2))
    for site_label,xyz in zip(site_labels,sites_cart):
        atom_entry=site_label.split('"')[1]
        name=atom_entry[0:4]
        resname=atom_entry[5:8]
        altloc=atom_entry[9]
        resnum=int(atom_entry[10:15])
        key = (resnum,resname,name)  
        assert key not in atom_xyz_dict, (atom_entry,key)
        atom_xyz_dict[key]=xyz
        #print(key,site_label,xyz)

    clashes_list=[]
    for f in filtered_dicts:
        # print(f.id_str())
        # print(f.as_selection_string())
        # # print(f.atom_selection)
        # # print(f.atoms_info)
        # # print(f.xyz)
        # # overlap = f.overlap - CLASH_SCAN_BUFFER
        # # vdw_sep= separation-overlap
        json_obj = json.loads(f.as_JSON())
        symop_str=json_obj["target_symop"]
        if symop_str is None:
            symop_str=""
        assert json_obj["symop"] is None,(type(json_obj["symop"]),json_obj["symop"])

        overlap = f.overlap + CLASH_SCAN_BUFFER
        
        #print(f.as_table_row_phenix()[:-1][1].split())

        site_labels = [] 
        for atom_entry in f.as_table_row_phenix()[:-1]:
            atom_entry=atom_entry[1:]
            altloc=atom_entry[0]
            resnum=int(atom_entry[1:6])
            resname=atom_entry[7:10]
            name=atom_entry[11:15]
            site_labels.append((resnum,resname,name))
        assert len(site_labels)==2
        coords= [atom_xyz_dict[key] for key in site_labels]
        vdw_sum= sep(*coords) - overlap
        #print(overlap)
        #print(sep(*coords))
        clashes_list.append((
            *site_labels,
            vdw_sum,
            symop_str
            )
        )
        assert vdw_sum >0,(sep(*coords),f.overlap,CLASH_SCAN_BUFFER,vdw_sum)
# 'pdb=" N  AARG A  62 "'
# '"pdb= A  62 AARG  C "'
    # model_statistics_geometry = model.geometry_statistics(
    #   use_hydrogens=None, condensed_probe=False, fast_clash=False)
    # model_statistics_geometry_result = \
    #   model_statistics_geometry.result()
    # ramalyze  = model_statistics_geometry_result.ramachandran.ramalyze
    # omegalyze = model_statistics_geometry_result.omega.omegalyze
    # rotalyze  = model_statistics_geometry_result.rotamer.rotalyze
    # cbetadev  = model_statistics_geometry_result.c_beta.cbetadev
    # clashes   = model_statistics_geometry_result.clash.clashes
    #print(clashes.clash_dict)
    #clashes.show()
    # probe_clashscore_manager(
    #     h_pdb_string=input_str,
    #     nuclear=nuclear,
    #     fast=self.fast,
    #     condensed_probe=self.condensed_probe,
    #     largest_occupancy=occ_max,
    #     b_factor_cutoff=b_factor_cutoff,
    #     use_segids=use_segids,
    #     verbose=verbose,
    #     model_id=model.id)


    #assert proxies_info_nonbonded is not None
    #nonbonded_list = proxies_info_nonbonded[0]



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
    ordered_atom_label_xyz_dict={}
    def sep(X,Y):
        return np.sqrt(np.sum((np.array(X)-np.array(Y))**2))
    for site_label,xyz in zip(og_site_labels,og_sites_cart):
        #key = site_label[:9]+site_label[10:]
        atom_entry=site_label.split('"')[1]
        name=atom_entry[0:4]
        resname=atom_entry[5:8]
        altloc=atom_entry[9]
        resnum=int(atom_entry[10:15])
        key = (resnum,resname,name)
        if key not in ordered_atom_label_xyz_dict:
            ordered_atom_label_xyz_dict[key]=[]
        ordered_atom_label_xyz_dict[key].append((site_label,xyz))
        #print(key,site_label,xyz)

####

    out_data=[]
    for i, item in enumerate(clashes_list):
        if i%100000==0:
            print(f"{i}/{len(clashes_list)}")
        site_label_A          = item[0]
        site_label_B          = item[1]
        vdw_sum        = item[2] # but for clashes
        symop_str      = item[3] 



        keyA=(new_resnum_to_old_resnum_dict[site_label_A[0]],)+site_label_A[1:]
        keyB=(new_resnum_to_old_resnum_dict[site_label_B[0]],)+site_label_B[1:]
        # if keyA=='pdb=" C  VAL A   1 "':
        #     if keyB[-4:-2]==" 2":
        #         print(keyA,keyB)
        #         print(site_labels[i_seq],site_labels[j_seq])
        distance_cutoff = 4
        
        debug_packing_added={}
        debug_nonpacking_added={}
        for (conformer_site_label_A,coordA) in ordered_atom_label_xyz_dict[keyA]:
            for (conformer_site_label_B,coordB) in ordered_atom_label_xyz_dict[keyB]:
                if sep(coordA,coordB) > distance_cutoff:
                    continue
                # if keyA=='pdb=" C  GLY A  59 "':
                #     if keyB[-5:-2]==" 60":
                #         print(model_distance,vdw_sum)
                #         print(conformer_site_label_A,conformer_site_label_B)
                #         print(site_labels[i_seq],site_labels[j_seq])
                crystal_packing_contact= "true" if symop_str.strip()!="" else "false"

                vdw_sum=round(vdw_sum,3)
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
                    tol=0.2
                    assert abs(debug_dict[debug_key]-vdw_sum)<tol, (datum,debug_dict[debug_key])
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