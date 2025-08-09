xyz_path=$1; 
out_tag=$2; shift 2

xyz_file=${xyz_path##*/}
xyz_handle=${xyz_file%.*}


out_handle=${xyz_handle}_${out_tag}

# defaults
hkl_name=refme
serial=999
wc=1
wu=1
wxc_scale=0.5
macro_cycles=5
calc_wE='false'
hold_water='false'
optimize_R='false'
shake=0

no_mlhl=true

while getopts ":o:u:c:n:s:d:w:h:r" flag; do
 case $flag in
    o) out_handle=$OPTARG
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
    ;;
    w) calc_wE='true'
    ;;
    h) hold_water='true'
    ;;
    r) optimize_R='true'
    ;;
   \?)
   echo INVALID FLAG
   ;;
 esac
done




echo $xyz_path $hkl_name $out_handle $wu $wc $macro_cycles $shake $calc_wE $hold_water $optimize_R


expected_path=$xyz_path
if [ ! -f $expected_path ]; then  # TODO checks after other files or make Refine and RptRefine do exit 0 on error 
    echo "File ${expected_path} not found!"
    exit 0
fi


# if ${shake}; then
#     paramFileTemplate=refine_water_hold_and_shake_protein_template.eff
# fi
paramFileTemplate=refine_water_bond_length_hold_template.eff
if $optimize_R; then 
  paramFileTemplate=refine_water_bond_length_hold_optimize_R_template.eff
fi
if $hold_water; then
  paramFileTemplate=refine_water_hold_template.eff
  if $optimize_R; then 
    paramFileTemplate=refine_water_hold_optimize_R_template.eff
  fi
fi
#paramFileTemplate=refine_water_bond_length_hold_template.eff
#paramFileTemplate=refine_water_hold_template_free_necessary_waters.eff

paramFile=${out_handle}_initial_refine.eff

cd $(dirname "$0")
rm tmp_refinement/$paramFile
cp $paramFileTemplate tmp_refinement/$paramFile
cp $xyz_path tmp_refinement/${xyz_handle}.pdb

cd tmp_refinement
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


logs_path="../../output/refine_logs"
if [ ! -d $logs_path ]; then
  mkdir $logs_path
fi

phenix.refine $paramFile > ../../output/refine_logs/${out_handle}.log


mv ${out_handle}_${serial}.pdb ../../output/${out_handle}.pdb  

cd .. 

if $calc_wE; then
  cd ../StructureGeneration
  bash GenerateHoltonData.sh ../output/${out_handle}.pdb  > HoltonScores/${out_handle}.log
fi








