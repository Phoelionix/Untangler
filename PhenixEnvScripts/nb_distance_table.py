
# from cctbx.geometry_restraints import nonbonded_distance_table
# from cctbx.geometry_restraints import nonbonded_deltas
# from cctbx import geometry_restraints

# p = geometry_restraints.nonbonded_params()
# d = p.distance_table

# d.setdefault("Si")["O"] = 1.5
# print(d["Si"]["O"])


# See phenix-2.0-5793/lib/python3.9/site-packages/cctbx/geometry_restraints/tst_process_nonbonded_proxies.py


# I don't know what I'm doing, just finding a place I can bootstrap from.

from cctbx.geometry_restraints import process_nonbonded_proxies
import mmtbx
import iotbx
import sys
from libtbx.utils import null_out

pdb_file = sys.argv[1]


params = mmtbx.model.manager.get_default_pdb_interpretation_params()
params.pdb_interpretation.allow_polymer_cross_special_position=True
params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
#pdb_inp = iotbx.pdb.input(lines=raw_records.split("\n"), source_info=None)
pdb_inp = iotbx.pdb.input(pdb_file)

model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out()
)

model.process(pdb_interpretation_params=params,
    make_restraints=True)


pnp_manager = process_nonbonded_proxies.manager(model=model)

grm = pnp_manager.model.get_restraints_manager().geometry
xrs = pnp_manager.model.get_xray_structure()
sites_cart  = pnp_manager.model.get_sites_cart()
site_labels = xrs.scatterers().extract_labels()
pair_proxies = grm.pair_proxies(
                    sites_cart  = sites_cart,
                    site_labels = site_labels)
proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted(
    by_value    = "delta",
    sites_cart  = sites_cart,
    site_labels = site_labels)


if proxies_info_nonbonded is not None:
    nonbonded_list = proxies_info_nonbonded[0]
else:
    assert False
    # nonbonded_list = []
    # # create 'empty' instance of results class
    # self._clashes = clashes(clashes_dict = dict())
    # self._hbonds  = hbonds(hbonds_dict = dict())
    # return

for item in nonbonded_list:
    i_seq          = item[1]
    j_seq          = item[2]
    model_distance = item[3]
    vdw_sum        = item[4]
    symop_str      = item[5] # TODO probably not necessary
    symop          = item[6]
    print(i_seq,j_seq,vdw_sum)

pair_proxies.nonbonded_proxies.show_histogram_of_model_distances(
    sites_cart=sites_cart,
    f=sys.stdout)

# That was nearly instant.
# Now need to just modify phenix to generate vdw sums for all possible altlocs or find the function that does

