
# from cctbx.geometry_restraints import nonbonded_distance_table
# from cctbx.geometry_restraints import nonbonded_deltas
# from cctbx import geometry_restraints

# p = geometry_restraints.nonbonded_params()
# d = p.distance_table

# d.setdefault("Si")["O"] = 1.5
# print(d["Si"]["O"])


# See phenix-2.0-5793/lib/python3.9/site-packages/cctbx/geometry_restraints/tst_process_nonbonded_proxies.py

from cctbx.geometry_restraints import process_nonbonded_proxies,pair_proxies
import mmtbx.model
import iotbx
import sys
from libtbx.utils import null_out
from scitbx_array_family_flex_ext import std_string
from scitbx.array_family import flex

pdb_file = sys.argv[1]


params = mmtbx.model.manager.get_default_pdb_interpretation_params()
params.pdb_interpretation.allow_polymer_cross_special_position=True
params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
#pdb_inp = iotbx.pdb.input(lines=raw_records.split("\n"), source_info=None)


tmp_pdb_file = "tmp_samealtloc.pdb"

with open(pdb_file) as f, open(tmp_pdb_file,"w") as w:
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
    model_input = iotbx.pdb.input(pdb_file),
    log         = null_out()
).get_xray_structure().scatterers().extract_labels()

num_skipped=0
for item in nonbonded_list:
    i_seq          = item[1]
    j_seq          = item[2]
    model_distance = item[3]
    vdw_sum        = item[4]
    symop_str      = item[5] # TODO probably not necessary
    symop          = item[6]

    # Skip if same atom but different altloc
    if og_site_labels[i_seq][:9]+og_site_labels[i_seq][10:]==og_site_labels[j_seq][:9]+og_site_labels[j_seq][10:]:
        num_skipped+=1    
        continue
    print(og_site_labels[i_seq],og_site_labels[j_seq],vdw_sum)
    #print(site_labels[i_seq],site_labels[j_seq],vdw_sum)
print(f"number of cross-conformer nonbonds: {len(nonbonded_list)-num_skipped}")

pair_proxies.nonbonded_proxies.show_histogram_of_model_distances(
    sites_cart=sites_cart,
    f=sys.stdout)

# That was nearly instant.
# Now need to just modify phenix to generate vdw sums for all possible altlocs or find the function that does

