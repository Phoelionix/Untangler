#%%

from __future__ import print_function
import sys 
import shutil
from iotbx.data_manager import DataManager
import mmtbx

def test_each_model(i_model,model_path):
  file_name = 'tst_vdw_and_h_bond_%02d.pdb' % (i_model+1)
  shutil.copy(model_path,file_name)
  dm = DataManager()
  dm.process_model_file(file_name)
  model = dm.get_model(file_name)

  useNeutronDistances = False
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = useNeutronDistances
  model.process(make_restraints=True)

  # Get information on VdW radii.
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib(use_neutron_distances = useNeutronDistances)
  ph = model.get_hierarchy()
  for atom in ph.atoms():
    ag = atom.parent()
    md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
      residue_name=ag.resname, atom_names=ag.atoms().extract_name())
    atom_dict = md.atom_dict()
    vdw = model.get_specific_vdw_radius(atom.i_seq)
    vdwh = model.get_specific_vdw_radius(atom.i_seq, vdw_radius_without_H=True)
    h_bond = model.get_specific_h_bond_type(atom.i_seq)
    ion_radius = model.get_specific_ion_radius(atom.i_seq)
    if ion_radius is None:
      ion_radius_str = '    -'
    else:
      ion_radius_str = '%7.2f' % ion_radius
    if vdwh is None:
      vdwh = -99
    # print('  %s  %7.2f %7.2f %s %s' % (atom.quote(),
    #                                    vdw,
    #                                    vdwh,
    #                                    h_bond,
    #                                    ion_radius_str))
    print('  %s  %7.2f %7.2f %s %s' % (atom.quote(),
                                       vdw,
                                       vdwh,
                                       h_bond,
                                       ion_radius_str))
    name = atom.name.strip()
    try:
      te = atom_dict[name].type_energy
    except:
      continue
    vdw_radius = ener_lib.lib_atom[te].vdw_radius
    assert vdw_radius==vdw
  print('OK')

def main():
  model_paths = sys.argv[1:]
  for i,f in enumerate(model_paths):
    test_each_model(i,f)

if __name__ == '__main__':
  main()

# %%
