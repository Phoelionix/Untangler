# bash Refine.sh ../data/4et8.pdb refined -d ../data/lys_galli_aug_real.mtz -n 5 -s 3 -gzt

# TODO folder with file containing commands used, shell file, and .eff file. 

xyz_path=$1; 
out_tag=$2; shift 2 # redundant?

xyz_file=${xyz_path##*/}
xyz_handle=${xyz_file%.*}


out_handle_override='false'

# defaults
hkl_name=refme
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
refine_alt_locs='false'
refine_no_hold='false'
no_mlhl=true
generate_r_free='false'
turn_off_bulk_solvent='false'

while getopts ":o:u:c:n:s:d:whpragtz" flag; do
 case $flag in
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
    d) hkl_name=$OPTARG
       hkl_name=${hkl_name##*/}
       hkl_name=${hkl_name%.*}
    ;;
    w) calc_wE='true'
    ;;
    h) hold_water='true'
    ;;
    p) hold_protein='true'
    ;;
    r) optimize_R='true'
    ;;
    a) refine_alt_locs='true'
    ;;
    g) generate_r_free='true'
    ;;
    t) turn_off_bulk_solvent='true'
    ;;
    z) refine_no_hold='true'
    ;;
   \?)
   echo INVALID FLAG
   ;;
 esac
done


if ! $out_handle_override; then
  out_handle=${xyz_handle}-${hkl_name}_${out_tag}

fi 

echo $xyz_path $hkl_name $out_handle $wu $wc $macro_cycles $shake $calc_wE $hold_water $optimize_R $refine_alt_locs $generate_r_free $refine_no_hold $turn_off_bulk_solvent


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

#paramFileTemplate=refine_water_bond_length_hold_template.eff
#paramFileTemplate=refine_water_hold_template_free_necessary_waters.eff

paramFile=${out_handle}_initial_refine.eff

cd $(dirname "$0")
mkdir -p tmp_refinement


mkdir -p tmp_refinement/$out_handle # Do refinement in own directory in attempt to stop seg faults when parallel. Possibly issue is due to the annoying .status.pkl file that is created

rm -f tmp_refinement/$out_handle/$paramFile
cp $paramFileTemplate tmp_refinement/$out_handle/$paramFile
cp $xyz_path tmp_refinement/$out_handle/${xyz_handle}.pdb


cd tmp_refinement/$out_handle


sed "s/XYZ_TEMPLATE/${xyz_handle}/g" $paramFile > tmp.$$
mv tmp.$$ $paramFile
sed  "s/HKL_TEMPLATE/${hkl_name}/g" $paramFile  > tmp.$$
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

if $no_mlhl; then
  sed "s/target = auto ml \*mlhl ml_sad ls mli/target = *auto ml mlhl ml_sad ls mli/g" $paramFile > tmp.$$ 
  mv tmp.$$ $paramFile
fi

if $refine_alt_locs; then 
  echo "Refining altlocs option not implemented (but easy to implement...)"
  exit
fi

if $generate_r_free; then
  sed  "s/generate = False/generate = True/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
fi

if $turn_off_bulk_solvent; then 
  sed "s/bulk_solvent_and_scale = True/bulk_solvent_and_scale = False/g" $paramFile  > tmp.$$ 
  mv tmp.$$ $paramFile
fi

logs_path="../../../output/refine_logs"
if [ ! -d $logs_path ]; then
  mkdir $logs_path
fi

# Broad sweep attempt to stop phenix segfaulting when run in parallel
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export TMPDIR=$PWD/tmp
mkdir -p $TMPDIR

phenix.refine $paramFile > $logs_path/${out_handle}.log
unset TMPDIR

mv ${out_handle}_${serial}.pdb ../../../output/${out_handle}.pdb  

cd ../.. 

if $calc_wE; then
  cd ../StructureGeneration
  bash GenerateHoltonData.sh ../output/${out_handle}.pdb  > HoltonScores/${out_handle}.log
fi








