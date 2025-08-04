#xyz_name=longrange_traps_refine_12_modelRptFreerWaters # In StructureGeneration/data
#out_handle=${xyz_name}_problemafix
xyz_name=$1; 
out_handle=${xyz_name}_refineSpecific
refinement_selection=$2; shift 2; # e.g. "resid 43 or resid 44 or resid 35"


hkl_name=refme # In StructureGeneration/data
serial=999
wc=0.01
wu=0.5
wxc_scale=0.5

macro_cycles=4

copy_to_data_folder='false'
while getopts ":o:u:c:n:d" flag; do
 case $flag in
    o) out_handle=$OPTARG
    ;;
    u) wu=$OPTARG
    ;;
    c) wc=$OPTARG
    ;;
    n) macro_cycles=$OPTARG
    ;;
    # s) shake='true'
    # ;;
    d) copy_to_data_folder='true'
    ;;
   \?)
   echo INVALID FLAG
   ;;
 esac
done

echo $refinement_selection $out_handle $macro_cycles $wu $wc $copy_to_data_folder

#paramFileTemplate=refine_water_hold_template.eff
paramFileTemplate=refine_only_problematic.eff
#paramFileTemplate=refine_water_hold_template_free_necessary_waters.eff

cd $(dirname "$0")
paramFile=${out_handle}_initial_refine.eff
rm tmp_refinement/$paramFile
cp $paramFileTemplate tmp_refinement/$paramFile

cd tmp_refinement
sed "s/XYZ_TEMPLATE/${xyz_name}/g" $paramFile > tmp.$$
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
######
sed "s/individual = (resid TEMPLATE or )/individual = ${refinement_selection}/g" $paramFile  > tmp.$$ 
mv tmp.$$ $paramFile

    

phenix.refine $paramFile


mv ${out_handle}_${serial}.pdb ../output/${out_handle}.pdb  

cd .. 

bash GenerateHoltonData.sh output/${out_handle}.pdb   > HoltonScores/${out_handle}.log


if $copy_to_data_folder; then
cp output/${out_handle}.pdb  ../data/${out_handle}.pdb
fi






