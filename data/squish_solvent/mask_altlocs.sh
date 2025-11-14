set -u 

# handle="8conf_waterSep"
# model="$handle.pdb"
# mtz="1ahoSF.mtz"

handle="8conf_waterSep_unrestrained0"
model="$handle.pdb"
mtz="$handle.mtz"
altlocs_to_keep="C"

#https://phenix-online.org/documentation/reference/reciprocal_space_arrays.html


rm -rf "out/"; mkdir -p out

mtz=$(realpath "$mtz")
out_path=our_altlocs_only.mtz

if [ -f $out_path ]; then 
    mv $out_path $out_path#
fi 

#cp $mtz out/recip.mtz
cp $model out/model.pdb

rm -f "out/recip.mtz"  # NB: If already exists, below command will get "multiple equally valid arrays found" and it will skip.
phenix.reciprocal_space_arrays $model $mtz output="out/recip.mtz"



out_path=$(realpath $out_path )


cd out
rm -f $out_path
subtract_altlocs.com $altlocs_to_keep model.pdb recip.mtz outfile=$out_path FP=FOBS SIGFP=SIGFOBS  FC_pair=FOBS,PHIFCALC


#altlocs=(`awk '/^ATOM|^HETAT/ {print substr($0,17,1)}' model.pdb | sort -u`)
#altlocs="${altlocs[@]}""
altlocs_array_string=$(echo "$altlocs_to_keep" | grep -o .)
echo $altlocs_array_string
awk -v altlocs_string="$altlocs_array_string" '

BEGIN {
    split(altlocs_string, altlocs_array, " ")
}

# if line starts with
/^(ATOM|HETATM)/ {
    conf=substr($0,17,1);
    for (i in altlocs_array)
        if (altlocs_array[i] == conf)
            print; next
}
# else
{ print }

' model.pdb > model_subset.pdb




bash ~/Untangler/Refinement/Refine.sh $(realpath model_subset.pdb) $out_path -n 0 -o mask_altlocs_test
bash ~/Untangler/Refinement/Refine.sh $(realpath model_subset.pdb) $mtz -n 0 -o mask_altlocs_control

# phenix.mtz.dump file
# phenix.mtz.dump --show_column_data file