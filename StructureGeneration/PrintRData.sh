
model=$1
pattern=$2

model_handle=$(basename "${1%.*}")
phenix.model_vs_data $model  $pattern #> $(dirname "$0")/output/SF_check_${model_handle}.log

