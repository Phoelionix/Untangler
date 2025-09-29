set -u 

parent_dir=../../Refinement/tmp_refinement/ 
handle=best_for_014_loopEnd15
dir=${parent_dir}${handle}/

cp ${dir}refmacout.pdb .
cp ${dir}refmacout.mtz .

add_waters_for_altlocs.com 

python3.9 delete_water_copies.py