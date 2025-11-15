# bash Refine.sh ../data/4et8.pdb refined -d ../data/lys_galli_aug_real.mtz -n 5 -s 3 -gzt

# TODO folder with file containing commands used, shell file, and .eff file. 

set -u

xyz_path=$1; 
hkl_path=$2; shift 2 # redundant?

xyz_file=${xyz_path##*/}
xyz_handle=${xyz_file%.*}
hkl_file=${hkl_path##*/}
hkl_handle=${hkl_file%.*}


out_handle_override='false'

# defaults
serial=999
wc=1
wu=1
wxc_scale=0.5
macro_cycles=5
calc_wE='false'
hold_water='false'
hold_protein='false'
optimize_R='false'
shake=0
refine_no_hold='false'
no_mlhl=true
generate_r_free='false'
turn_off_bulk_solvent='false'
restrain_movement='true'
refine_occupancies='false'
ordered_solvent='false'
disable_ADP='false'
disable_CDL='false' # Disable conformation-dependent library
altlocs_to_refine=''

max_sigma_movement_water=0.1

refine_hydrogens='false'

while getopts ":a:o:u:c:n:s:q:whprgtzACHORS" flag; do
 case $flag in
    a) altlocs_to_refine=$OPTARG
    ;;
    o) out_handle=$OPTARG
       out_handle_override='true'
    ;;
    u) wu=$OPTARG
    ;;
    c) wc=$OPTARG
    ;;
    n) macro_cycles=$OPTARG
    ;;
    s) shake=$OPTARG
    ;;
    q) max_sigma_movement_water=$OPTARG
    ;;
    w) calc_wE='true'
    ;;
    h) hold_water='true'
    ;;
    p) hold_protein='true'
    ;;
    r) optimize_R='true'
    ;;
    g) generate_r_free='true'
    ;;
    t) turn_off_bulk_solvent='true'
    ;;
    z) refine_no_hold='true'
    ;;
    A) disable_ADP='true'
    ;;
    C) disable_CDL='true'
    ;;
    H) refine_hydrogens='true'
    ;;
    O) refine_occupancies='true'
    ;;
    R) restrain_movement='false' 
    ;;
    S) ordered_solvent='true'
    ;;
   \?)
   echo INVALID FLAG
   ;;
 esac
done


if ! $out_handle_override; then
  out_handle=${xyz_handle}-${hkl_handle}
  if [ $altlocs_to_refine != "" ]; then
    out_handle=${out_handle}-${altlocs_to_refine}
  fi
fi 

echo $xyz_path $hkl_path $out_handle $wu $wc $macro_cycles $shake $calc_wE $hold_water $optimize_R $generate_r_free $refine_no_hold $turn_off_bulk_solvent $restrain_movement $refine_hydrogens $refine_occupancies 


expected_path=$xyz_path
if [ ! -f $expected_path ]; then  # TODO checks after other files or make Refine and RptRefine do exit 0 on error 
    echo "File ${expected_path} not found!"
    exit 0
fi


# if ${shake}; then
#     paramFileTemplate=refine_water_hold_and_shake_protein_template.eff
# fi


#paramFileTemplate=refine_water_bond_length_hold_template.eff
paramFileTemplate=refine_no_hold_template.eff
#paramFileTemplate=refine_hold_altlocs.eff
if $optimize_R; then 
  #paramFileTemplate=refine_water_bond_length_hold_optimize_R_template.eff
  paramFileTemplate=refine_no_hold_optimize_R_template.eff
fi
if $hold_water; then
  paramFileTemplate=refine_water_hold_template.eff
  if $optimize_R; then 
    paramFileTemplate=refine_water_hold_optimize_R_template.eff
  fi
fi

if $hold_protein; then
  paramFileTemplate=refine_protein_hold_template.eff
fi

#TEMPORARY
if $refine_no_hold; then
  paramFileTemplate=refine_no_hold_template.eff
fi

echo  $paramFileTemplate

#paramFileTemplate=refine_water_bond_length_hold_template.eff
#paramFileTemplate=refine_water_hold_template_free_necessary_waters.eff

paramFile=${out_handle}_in.eff



#xyz_path=$(realpath -s --relative-to="$(dirname "$0")" "$xyz_path" )
#hkl_path=$(realpath -s --relative-to="$(dirname "$0")" "$hkl_path" )
xyz_path=$(realpath "$xyz_path" )
hkl_path=$(realpath "$hkl_path" )

cd $(dirname "$0")
mkdir -p tmp_refinement



mkdir -p tmp_refinement/$out_handle # Do refinement in own directory in attempt to stop seg faults when parallel. Possibly issue is due to the annoying .status.pkl file that is created

rm -f tmp_refinement/$out_handle/$paramFile
cp $paramFileTemplate tmp_refinement/$out_handle/$paramFile
cp $xyz_path tmp_refinement/$out_handle/${xyz_handle}.pdb


cd tmp_refinement/$out_handle



xyz_subset_handle=$xyz_handle
if [ $altlocs_to_refine != "" ]; then
  xyz_subset_handle=${xyz_handle}-${altlocs_to_refine}

  echo "Masking out other altlocs"
  # (Approximately) remove contributions of other altlocs from Fobs. 
  new_hkl_path=$(realpath "altlocs_masked.mtz")
  tmp_dir=$(realpath "mask_out/")
  rm -rf $tmp_dir; mkdir -p $tmp_dir
  bash $(dirname "$0")/mask_altlocs.sh ${xyz_handle}.pdb $hkl_path $altlocs_to_refine $new_hkl_path $tmp_dir &> /dev/null
  # Get structure file of the atoms with the specified altloc labels
  bash $(dirname "$0")/make_altloc_subset.sh ${xyz_handle}.pdb $altlocs_to_refine ${xyz_subset_handle}.pdb &> /dev/null
 
fi



sed "s/XYZ_TEMPLATE/${xyz_subset_handle}/g" $paramFile > tmp.$$
mv tmp.$$ $paramFile
sed  "s@HKL_TEMPLATE_PATH@${hkl_path}@g" $paramFile  > tmp.$$
mv tmp.$$ $paramFile
sed  "s/PREFIX_TEMPLATE/${out_handle}/g" $paramFile  > tmp.$$
mv tmp.$$ $paramFile
sed  "s/serial = None/serial = ${serial}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile
sed  "s/serial = None/serial = ${serial}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile
sed  "s/    wc = 1/    wc = ${wc}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile
sed  "s/    wu = 1/    wu = ${wu}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile
sed  "s/NUM_MACRO_CYCLES/${macro_cycles}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile
sed  "s/wxc_scale = 0.5/wxc_scale = ${wxc_scale}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile
sed  "s/SHAKE_TEMPLATE/${shake}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile
sed  "s/      sigma = 0.1/      sigma = ${max_sigma_movement_water}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile

if $no_mlhl; then
  sed "s/target = auto ml \*mlhl ml_sad ls mli/target = *auto ml mlhl ml_sad ls mli/g" $paramFile > tmp.$$ 
  mv tmp.$$ $paramFile
fi

if $refine_occupancies; then 
  sed  "s/tls occupancies/tls *occupancies/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
  sed  "s/remove_selection = All/remove_selection = None/g" $paramFile  > tmp.$$
  mv tmp.$$ $paramFile
fi

if $generate_r_free; then
  sed  "s/generate = False/generate = True/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
fi

if $turn_off_bulk_solvent; then 
  sed "s/bulk_solvent_and_scale = True/bulk_solvent_and_scale = False/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
fi

if $refine_hydrogens; then 
  sed "s/refine = individual \*riding Auto/refine = *individual riding Auto/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
  sed "s/xh_bond_distance_deviation_limit = 0/xh_bond_distance_deviation_limit = 999/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
  # sed "s/real_space_optimize_x_h_orientation = True/real_space_optimize_x_h_orientation = False/g" $paramFile  > tmp.$$ 
  # mv tmp.$$ $paramFile
fi

if $ordered_solvent; then 
  sed "s/ordered_solvent = False/ordered_solvent = True/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
fi

if $disable_ADP; then 
  sed "s/\*individual_adp/individual_adp/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
fi

if $disable_CDL; then 
  sed "s/cdl = True/cdl = False/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
fi

logs_path="../../../output/refine_logs"
mkdir -p $logs_path

if ! $restrain_movement; then 
    sed 's/reference_coordinate_restraints {\n      enabled = True/reference_coordinate_restraints {\n      enabled = False/g' $paramFile  > tmp.$$ 
    mv tmp.$$ $paramFile
fi 




# Broad sweep attempt to stop phenix segfaulting when run in parallel
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export TMPDIR=$PWD/tmp_${out_handle}
mkdir -p $TMPDIR

# env -i PATH=/usr/local/phenix-2/build/bin:/usr/bin:/bin \
#   PHENIX=/usr/local/phenix-2 \
echo "Refining $out_handle"
phenix.refine $paramFile > $logs_path/${out_handle}.log
unset TMPDIR

final_structure=${out_handle}_${serial}.pdb

if [[ $altlocs_to_refine != "" ]]; then
  # Add back in the conformations we didn't refine.
  final_structure=${xyz_handle}_updated.pdb  
  bash $(dirname "$0")/update_altloc_subset.sh $xyz_handle.pdb ${out_handle}_${serial}.pdb  $final_structure
fi

out_path=$(realpath ../../../output/${out_handle}.pdb  )

cp $final_structure $out_path 

cd ../.. 

if $calc_wE; then
  echo "Calculating wE for $out_path"
  cd ../StructureGeneration
  echo $(realpath -s --relative-to="./" "$out_path" )
  bash GenerateHoltonData.sh $(realpath -s --relative-to="./" "$out_path" )  > HoltonScores/${out_handle}.log
fi








