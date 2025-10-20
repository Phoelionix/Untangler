set -u 

handle="5conf_loopEnd4-4"
model="$handle.pdb"
mtz="$handle.mtz"

#https://phenix-online.org/documentation/reference/reciprocal_space_arrays.html

mkdir -p out

if [ -f refme_minusol.mtz ]; then 
    mv refme_minusol.mtz refme_minusol#.mtz
fi 

rm -f "out/recip.mtz"
phenix.reciprocal_space_arrays $model $mtz output="out/recip.mtz"

phenix.reduce $model -trim > out/noH.pdb
cp $model out/model.pdb

cd out
squish_solvent_runme.com model.pdb recip.mtz FP=FOBS SIGFP=SIGFOBS FC_subtract=FMODEL,PHIFMODEL FC_addback=FCALC,PHIFCALC
#squish_solvent_runme.com model.pdb recip.mtz FC_subtract=FOBS,PHIFMODEL FC_addback=FOBS,PHIFCALC

mv refme_minusol.mtz ../
