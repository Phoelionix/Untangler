set -u 

# handle="8conf_waterSep"
# model="$handle.pdb"
# mtz="1ahoSF.mtz"


# READ THIS
# Final model contains:
# 1. all atoms in $subset_model
# 2. all atoms in $model **whose altlocs are NOT in $subset_model**

model=$1
subset_model=$2
out_path=$3



model_altlocs=`
cat $model |\
awk '! /^ATOM|^HETAT/{next}\
 {conf=substr($0,17,1);occ=substr($0,55,6);\
  if(conf==" ")next;\
  ++count[conf];sum[conf]+=occ}\
 ! order[conf]{++c;order[conf]=c}\
 END{for(conf in count){\
   print order[conf],conf,sum[conf]/count[conf],count[conf]}}' |\
sort -g |\
while read -r cnfnum cnfchar rest; do
    echo $cnfchar
done
`

subset_altlocs=`
cat $subset_model |\
awk '! /^ATOM|^HETAT/{next}\
 {conf=substr($0,17,1);occ=substr($0,55,6);\
  if(conf==" ")next;\
  ++count[conf];sum[conf]+=occ}\
 ! order[conf]{++c;order[conf]=c}\
 END{for(conf in count){\
   print order[conf],conf,sum[conf]/count[conf],count[conf]}}' |\
sort -g |\
while read -r cnfnum cnfchar rest; do
    #subset_altlocs="$subset_altlocs$test "
    echo $cnfchar
done
`

subset_altlocs=$(echo $subset_altlocs | sed 's/ //g')
model_altlocs_kept=""
for m in $model_altlocs; do 
    if [[ ! $m =~ [$subset_altlocs] ]]; then
        model_altlocs_kept=${model_altlocs_kept}${m}
    fi
done



python3.9 $(dirname "$0")/../StructureGeneration/CombineUpdateEnsembles.py $out_path $model $model_altlocs_kept $subset_model $subset_altlocs



# awk -v model_altlocs_string="$model_altlocs" exclude_altlocs_string="$subset_altlocs" '

# BEGIN {
#     split(model_altlocs_string, model_altlocs_array, " ")
#     split(exclude_altlocs_string, exclude_altlocs_array, " ")
# }

# # if line starts with
# /^(ATOM|HETATM)/ {
#     conf=substr($0,17,1);
#     for (i in altlocs_array);
#         if (altlocs_array[i] == conf)
# }
# # else
# { print }

# ' 'FNR==NR { a[FNR""] = $0; next } { print a[FNR""], $0 }' $model $subset_model > $out_path



# phenix.mtz.dump file
# phenix.mtz.dump --show_column_data file