set -u

# model=$1; 
# data=$2
#model=/home/speno/Untangler/output/lrt_shaken_reconnected.pdb
#model=/home/speno/Untangler/output/cheat_lrt_shaken.pdb
#model=/home/speno/Untangler/output/lrt_shaken_connected_to_match_huge_shake.pdb
#model=/home/speno/Untangler/output/lrt_mixed_up.pdb
#model=/home/speno/Untangler/output/lrt_mixed_up_randomized.pdb
model=/home/speno/Untangler/output/lrt_mixed_up_reconnected.pdb
#model=/home/speno/Untangler/output/cheat_lrt_mixed_up.pdb
#model=/home/speno/Untangler/output/lrt_mixed_up_Cys12Swap.pdb
data=/home/speno/Untangler/data/synthetic/refme.mtz;
#model=/home/speno/Untangler/output/44PS_3conf_mixedup.pdb
#model=/home/speno/Untangler/output/cheat_44PS_3conf_mixedup.pdb
#data=/home/speno/Untangler/data/synthetic/44PS_3conf.mtz

xyz_file=${model##*/}
xyz_handle=${xyz_file%.*}
hkl_file=${data##*/}
hkl_handle=${hkl_file%.*}




# Higher weights at first mimic phenix.refine behavior
bash /home/speno/Untangler/Refinement/Refine.sh $model $data -c 12 -u 50 -n 1 -q 0.1 -H
bash /home/speno/Untangler/Refinement/Refine.sh output/$xyz_handle-${hkl_handle}.pdb $data -c 9 -u 40 -n 1 -q 0.1 -H
for i in $(seq 2 15); do 
    bash /home/speno/Untangler/Refinement/Refine.sh output/$xyz_handle-${hkl_handle}$i.pdb $data -c 1 -u 1 -n 1 -q 0.1 -H
done


