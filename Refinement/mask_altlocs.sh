set -u 

# handle="8conf_waterSep"
# model="$handle.pdb"
# mtz="1ahoSF.mtz"


model=$1
mtz=$2
altlocs_to_keep=$3
out_path=$4
working_dir=$5

#https://phenix-online.org/documentation/reference/reciprocal_space_arrays.html


if [ -f $out_path ]; then 
    mv $out_path $out_path#
fi 

cp $model $working_dir/model.pdb

rm -f "$working_dir/recip.mtz"  # NB: If already exists, below command will get "multiple equally valid arrays found" and it will skip.
phenix.reciprocal_space_arrays $model $mtz output="$working_dir/recip.mtz"



out_path=$(realpath $out_path )

cd $working_dir
rm -f $out_path
subtract_altlocs.com $altlocs_to_keep model.pdb recip.mtz outfile=$out_path FP=FOBS SIGFP=SIGFOBS  FC_pair=FOBS,PHIFCALC