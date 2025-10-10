#! /bin/tcsh -f
set -u 

handle="refmacout"
model="$handle.pdb"
mtz="$handle.mtz"

#https://phenix-online.org/documentation/reference/reciprocal_space_arrays.html

mkdir -p out

if [ -f refme_minusol.mtz ]; then 
    mv refme_minusol.mtz refme_minusol#.mtz
fi 

phenix.reduce $model -trim > out/noH.pdb
cp noH.pdb out/
cp $mtz out/refmacout.mtz

cd out
squish_solvent_runme.com

mv refme_minusol.mtz ../

