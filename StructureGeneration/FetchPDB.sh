
pdbID=$1



wget -O $pdbID.cif https://files.rcsb.org/download/$pdbID.cif

wget -O $pdbID.pdb https://files.rcsb.org/download/$pdbID.pdb

wget -O $pdbID-sf.cif https://files.rcsb.org/download/$pdbID-sf.cif

phenix.cif_as_mtz $pdbID-sf.cif