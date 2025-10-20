from Bio.PDB import PDBParser,Structure
from Bio.PDB.Atom import Atom,DisorderedAtom
from Bio.PDB.Residue import Residue # Note we don't want the disorderedresidue here, it refers to different residue types for same res seq num.
from Bio.PDB.StructureBuilder import StructureBuilder
from pulp import *
import UntangleFunctions 
from LinearOptimizer.Swapper import Swapper
import os
import random
import numpy as np
import itertools
from multiprocessing import Pool
import scipy.stats as st
from statistics import NormalDist
from LinearOptimizer.VariableID import *

header='''attribute: correlationCoefficient
match mode: 1-to-1
recipient: atoms'''

atoms:list[Atom] = PDBParser(reference_model).get_structure().get_atoms()

for atom in atoms:
	snum =
	:1	0.8280