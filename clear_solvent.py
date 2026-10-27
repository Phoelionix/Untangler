from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from UntangleFunctions import prepare_pdb
import sys

assert len(sys.argv)==3
pdb,out_path=sys.argv[1:]

prepare_pdb(pdb,out_path,allow_no_altloc=True,keep_anisou=True)
ordered_atom_lookup=OrderedAtomLookup(out_path,waters=True)

ordered_atom_lookup.remove_water_around_sidechains(radius=1.5,same_altlocs=False)
ordered_atom_lookup.output_as_pdb_file(out_path,out_path)