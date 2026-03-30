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
random_seed=42
serial=999
wc=1
wu=1
max_ordered_solvent_B=80
const_shrink_donor_acceptor=0
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
disable_movement_restraint='false'
refine_occupancies='false'
refine_water_occupancies='false'
ordered_solvent='false'
disable_ADP='false'
group_adp='false'
ADP_only='false'
water_and_H_only='false'
disable_CDL='false' # Disable conformation-dependent library
disable_nqh_flips='false'
altlocs_to_refine=''
user_param_file=''
filter_ordered_solvent='false'
clear_out_solvent_mode='false'
real_space_refine='false'
hold_main_chain='false'

max_sigma_movement_restraint=0.1

refine_hydrogens='false'

restrain_movement_of_protein='false' # Note this does nothing when True if disable_movement_restraint is True

fixed_water_occupancy='false' # Fix water occupancies at value of ordered_solvent_occupancy

reTry_on_fail='false' # You should not have any need to use this.


while getopts ":a:f:o:u:c:e:n:s:q:whprgtzACDFGHLMNOPRSTWXZ" flag; do
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
    e) user_param_file=$OPTARG
    ;;
    f) ordered_solvent_occupancy=$OPTARG
       fixed_water_occupancy='true' 
    ;;
    n) macro_cycles=$OPTARG
    ;;
    s) shake=$OPTARG
    ;;
    q) max_sigma_movement_restraint=$OPTARG
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
    D) ADP_only='true'
    ;;
    F) filter_ordered_solvent='true'
    ;;
    G) group_adp='true'
    ;;
    H) refine_hydrogens='true'
    ;;
    L) clear_out_solvent_mode='true'
    ;;
    M) hold_main_chain='true'
    ;;
    N) disable_nqh_flips='true'
    ;;
    O) refine_occupancies='true'
    ;;
    P) restrain_movement_of_protein='true'
    ;;
    R) disable_movement_restraint='true' 
    ;;
    S) ordered_solvent='true'
    ;;
    T) reTry_on_fail='true'
    ;;
    W) refine_water_occupancies='true'
    ;;
    X) real_space_refine='true'
    ;;
    Z) water_and_H_only='true'
    ;;
   \?)
   echo INVALID FLAG $flag
   ;;
 esac
done


if ! $out_handle_override; then
  out_handle=${xyz_handle}-${hkl_handle}
  if [ -n "$altlocs_to_refine" ]; then
    out_handle=${out_handle}-${altlocs_to_refine}
  fi
fi 

#echo $xyz_path $hkl_path $out_handle $wu $wc $macro_cycles $shake $calc_wE $hold_water $optimize_R $generate_r_free $refine_no_hold $turn_off_bulk_solvent $disable_movement_restraint $refine_hydrogens $refine_occupancies 


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
if [ -n "$altlocs_to_refine" ]; then
  xyz_subset_handle=${xyz_handle}-${altlocs_to_refine}

  echo "Masking out other altlocs"
  # (Approximately) remove contributions of other altlocs from Fobs. 
  new_hkl_path=$(realpath "altlocs_masked.mtz")
  tmp_dir=$(realpath "mask_out/")
  rm -rf $tmp_dir; mkdir -p $tmp_dir
  bash $(dirname "$0")/mask_altlocs.sh ${xyz_handle}.pdb $hkl_path $altlocs_to_refine $new_hkl_path $tmp_dir &> /dev/null
  # Get structure file of the atoms with the specified altloc labels
  bash $(dirname "$0")/../StructureGeneration/make_altloc_subset.sh ${xyz_handle}.pdb $altlocs_to_refine ${xyz_subset_handle}.pdb > /dev/null
  hkl_path=$new_hkl_path
 
fi

export TMPDIR=$PWD/tmp_${out_handle}
mkdir -p $TMPDIR

tmpfile=$TMPDIR/tmp.$$

sed "s/XYZ_TEMPLATE/${xyz_subset_handle}/g" $paramFile > $tmpfile
mv $tmpfile $paramFile
sed  "s@HKL_TEMPLATE_PATH@${hkl_path}@g" $paramFile  > $tmpfile
mv $tmpfile $paramFile
sed  "s/PREFIX_TEMPLATE/${out_handle}/g" $paramFile  > $tmpfile
mv $tmpfile $paramFile
sed  "s/serial = None/serial = ${serial}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
sed  "s/serial = None/serial = ${serial}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
sed  "s/    wc = 1/    wc = ${wc}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
sed  "s/    wu = 1/    wu = ${wu}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
sed  "s/NUM_MACRO_CYCLES/${macro_cycles}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
sed  "s/wxc_scale = 0.5/wxc_scale = ${wxc_scale}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
sed  "s/SHAKE_TEMPLATE/${shake}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
sed  "s/      sigma = 0.1/      sigma = ${max_sigma_movement_restraint}/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile


if $refine_hydrogens; then
  if $hold_protein; then 
    #echo "Hydrogens are the only protein atoms that will be refined"
    sed "s/individual = water/individual = water or element H/g" $paramFile > $tmpfile 
    mv $tmpfile $paramFile
  fi
fi

if $water_and_H_only; then 
    sed "s/individual = TEMPLATE_SITES_INDIVIDUAL/individual = water or element H/g" $paramFile > $tmpfile 
    mv $tmpfile $paramFile
fi

if $hold_main_chain; then 
    sed "s/individual = TEMPLATE_SITES_INDIVIDUAL/individual = not (name C or name CA or name N or name N1 or name NH or name O)/g" $paramFile > $tmpfile 
    mv $tmpfile $paramFile
fi

if $no_mlhl; then
  sed "s/target = auto ml \*mlhl ml_sad ls mli/target = *auto ml mlhl ml_sad ls mli/g" $paramFile > $tmpfile 
  mv $tmpfile $paramFile
fi

if $refine_occupancies; then 
  sed  "s/tls occupancies/tls *occupancies/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed  "s/remove_selection = All/remove_selection = None/g" $paramFile  > $tmpfile
  mv $tmpfile $paramFile
fi

if $refine_water_occupancies; then 
  sed  "s/tls occupancies/tls *occupancies/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed  "s/remove_selection = All/remove_selection = not water/g" $paramFile  > $tmpfile
  mv $tmpfile $paramFile
fi

if $generate_r_free; then
  sed  "s/generate = False/generate = True/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

if $turn_off_bulk_solvent; then 
  sed "s/bulk_solvent_and_scale = True/bulk_solvent_and_scale = False/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

if $refine_hydrogens; then 
  sed "s/refine = individual \*riding Auto/refine = *individual riding Auto/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed "s/xh_bond_distance_deviation_limit = 0/xh_bond_distance_deviation_limit = 999/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  # sed "s/real_space_optimize_x_h_orientation = True/real_space_optimize_x_h_orientation = False/g" $paramFile  > $tmpfile 
  # mv $tmpfile $paramFile
fi

if $real_space_refine; then 
  sed "s/individual_sites_real_space/*individual_sites_real_space/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

if $ordered_solvent; then 
  sed "s/ordered_solvent = False/ordered_solvent = True/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

if $group_adp; then 
  sed "s/\*individual_adp group_adp/individual_adp *group_adp/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi 

if $ADP_only; then 
  if $disable_ADP; then
    echo "ADP only but ADP disabled!"
    exit 
  fi 
  sed "s/\*individual_sites/individual_sites/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi 

if $disable_ADP; then 
  sed "s/\*individual_adp/individual_adp/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

if $disable_CDL; then 
  sed "s/cdl = True/cdl = False/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

logs_path="../../../output/refine_logs"
mkdir -p $logs_path

if $disable_movement_restraint; then 
  sed 's/reference_coordinate_restraints {\n      enabled = True/reference_coordinate_restraints {\n      enabled = False/g' $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi 

if $restrain_movement_of_protein; then 
  sed "s/selection = water/selection = all/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

if $disable_nqh_flips; then 
  sed "s/nqh_flips = True/nqh_flips = False/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

if $fixed_water_occupancy; then  
  sed "s/occupancy_min = 0.02/occupancy_min = $ordered_solvent_occupancy/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed "s/occupancy_max = 1.0/occupancy_max = $ordered_solvent_occupancy/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed "s/occupancy = 0.33/occupancy = $ordered_solvent_occupancy/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi

sed "s/const_shrink_donor_acceptor = 0/const_shrink_donor_acceptor = $const_shrink_donor_acceptor/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile


sed "s/individual = TEMPLATE_SITES_INDIVIDUAL/individual = None/g" $paramFile > $tmpfile 
mv $tmpfile $paramFile


if $filter_ordered_solvent; then
  sed 's/mode = \*second_half filter_only every_macro_cycle every_macro_cycle_after_first/mode = second_half *filter_only every_macro_cycle every_macro_cycle_after_first/g' $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed "s/ordered_solvent = False/ordered_solvent = True/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi
if $clear_out_solvent_mode; then
  sed 's/mode = \*second_half filter_only every_macro_cycle every_macro_cycle_after_first/mode = second_half *filter_only every_macro_cycle every_macro_cycle_after_first/g' $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed "s/ordered_solvent = False/ordered_solvent = True/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed "s/dist_min = 1.8/dist_min = 1.8/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
  sed "s/dist_min_altloc = 0.03/dist_min_altloc = 1.8/g" $paramFile  > $tmpfile 
  mv $tmpfile $paramFile
fi
sed "s/b_iso_max = 80.0/b_iso_max = $max_ordered_solvent_B/g" $paramFile  > $tmpfile 
mv $tmpfile $paramFile
# Broad sweep attempt to stop phenix segfaulting when run in parallel
# export OMP_NUM_THREADS=1
# export OPENBLAS_NUM_THREADS=1
# export MKL_NUM_THREADS=1
# export VECLIB_MAXIMUM_THREADS=1
# export NUMEXPR_NUM_THREADS=1

# export PYTHONFAULTHANDLER=1

# env -i PATH=/usr/local/phenix-2/build/bin:/usr/bin:/bin \
#   PHENIX=/usr/local/phenix-2 \
final_structure=${out_handle}_${serial}.pdb
if [ -f $final_structure ]; then 
  mv $final_structure $final_structure#
fi
echo "Refining $out_handle"
#echo $(realpath $paramFile) 
export PYTHONNOUSERSITE=1

#phenix.refine $paramFile > $logs_path/${out_handle}.log
error_file=$logs_path/${out_handle}_err.log
log_file=$logs_path/${out_handle}.log


num_attempts=1
if $reTry_on_fail; then 
  num_attempts=20 
fi

i=0
while true; do
  #user_param_file=/home/speno/Untangler/ConformationTree/output/split_conformations_restraints.eff
  #user_param_file=/home/speno/Untangler/ConformationTree/output/conformation_tree_restraints-4PSS_2conf6conf.eff # TEMPORARY
  #user_param_file=/home/speno/Untangler/ConformationTree/output/split_conformations_restraints-4PSS_2conf6conf_noWater.eff 
  #user_param_file=/home/speno/Untangler/ConformationTree/output/split_conformations_restraints-cov63_2conf6conf.eff
  #user_param_file=/home/speno/Untangler/ConformationTree/output/split_conformations_restraints-cov63_2conf6conf_noWater.eff
  #user_param_file=/home/speno/Untangler/ConformationTree/output/split_conformations_restraints_noNB-4PSS_6conf18conf.eff
  #user_param_file=/home/speno/Untangler/ConformationTree/output/split_conformations_restraints_noNB-4PSS_6conf18conf_reduced.eff
  #user_param_file=/home/speno/Untangler/Refinement/4PSS_group_adp.eff

  # XXX Assuming user param file is only ever conformation tree restraints...
  if (( $wc == 0 )); then  
    user_param_file=""
  fi


  failed=false
  phenix.refine main.random_seed=$random_seed $paramFile $user_param_file  2>$error_file 1> $log_file

  if [ -s $error_file ]; then
    failed=true
  fi

  if grep -q "Traceback (most recent call last)" "$log_file"; then
    failed=true
  fi
  if [ ! -f $final_structure ]; then 
    failed=true
    echo "$final_structure missing, refinement failed for unknown reason"
  fi
  if $failed; then
    mv $error_file $error_file#
    echo $error_file#
    i=$((i+1))
    echo
    if (( i >= num_attempts )); then 
      echo "Failed"
      exit 1 
    fi
    echo "Failed, retrying..."
  else
    rm $error_file
    break
  fi 
done
unset TMPDIR


if [ ! -f $final_structure ]; then 
  echo "$final_structure missing, refinement failed for unknown reason"
  exit
fi


if [ -n "$altlocs_to_refine" ]; then
  # Add back in the conformations we didn't refine.
  final_structure=${xyz_handle}_updated.pdb  
  bash $(dirname "$0")/update_altloc_subset.sh $xyz_handle.pdb ${out_handle}_${serial}.pdb  $final_structure
fi

out_path=$(realpath ../../../output/${out_handle}.pdb  )

cp $final_structure $out_path; cp ${out_handle}_${serial}.mtz  $(realpath ../../../output/${out_handle}.mtz  )

cd ../.. 

if $calc_wE; then
  echo "Calculating wE for $out_path"
  cd ../StructureGeneration
  bash GenerateHoltonData.sh $(realpath -s --relative-to="./" "$out_path" )  > HoltonScores/${out_handle}.log
fi








