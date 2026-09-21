set -u

# model=$1; 
# data=$2
model=/home/speno/Untangler/output/cheat_longrangetraps_start.pdb
data=/home/speno/Untangler/data/synthetic/refme.mtz;

xyz_file=${model##*/}
xyz_handle=${xyz_file%.*}
hkl_file=${data##*/}
hkl_handle=${hkl_file%.*}




# Higher weights at first mimic phenix.refine behavior
bash /home/speno/Untangler/Refinement/Refine.sh $model $data -c 12 -u 50 -n 1 -q 0.1 -H
bash /home/speno/Untangler/Refinement/Refine.sh output/$xyz_handle-${hkl_handle}.pdb $data -c 9 -u 40 -n 1 -q 0.1 -H
for i in $(seq 2 12); do 
    bash /home/speno/Untangler/Refinement/Refine.sh output/$xyz_handle-${hkl_handle}$i.pdb $data -c 1 -u 1 -n 1 -q 0.1 -H
done


