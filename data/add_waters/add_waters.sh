set -u 

rm -f new.pdb

mode="PHENIX"

parent_dir=../../Refinement/tmp_refinement/ 
#handle=best_for_014_loopEnd15
#handle=best_for_014_WaddedMod_loopEnd1
#handle=12conf_loopEnd4
#handle=4conf_refined
#handle=4conf_ref1_refined
handle=7conf_loopEnd4-4
dir=${parent_dir}${handle}/

if [ $mode == "REFMAC" ]; then

    cp ${dir}refmacout.pdb .
    cp ${dir}refmacout.mtz .

    add_waters_for_altlocs.com 
fi
if [ $mode == "PHENIX" ]; then 
    cp ${dir}${handle}_999.pdb .
    cp ${dir}${handle}_999.mtz .
    add_waters_for_altlocs.com ${dir}${handle}_999.pdb ${dir}${handle}_999.mtz
fi


python3.9 delete_water_copies.py new.pdb new.pdb
#python3.9 separate_altloc_set.py new.pdb new.pdb
#python3.9 separate_clashing_waters.py new.pdb new.pdb 1.8