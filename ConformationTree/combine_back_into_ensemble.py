import sys,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from LinearOptimizer.Tag import *
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler
from LinearOptimizer.Input import LP_Input
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from Bio.PDB import PDBParser,Structure,PDBIO
from UntangleFunctions import parse_symmetries_from_pdb



# Create full conformations from each altloc by copying the parent altlocs.