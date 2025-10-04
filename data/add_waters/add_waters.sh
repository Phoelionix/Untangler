set -u 

rm -f new.pdb

parent_dir=../../Refinement/tmp_refinement/ 
#handle=best_for_014_loopEnd15
#handle=best_for_014_WaddedMod_loopEnd1
handle=best_for_014_WaddedMod_initial
dir=${parent_dir}${handle}/

cp ${dir}refmacout.pdb .
cp ${dir}refmacout.mtz .

add_waters_for_altlocs.com 

python3.9 delete_water_copies.py new.pdb new.pdb