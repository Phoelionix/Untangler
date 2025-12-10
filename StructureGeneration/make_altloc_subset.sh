set -u 

# handle="8conf_waterSep"
# model="$handle.pdb"
# mtz="1ahoSF.mtz"


model=$1
altlocs_to_keep=$2
out_path=$3


altlocs_array_string=$(echo "$altlocs_to_keep" | grep -o .)

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

' $model > $out_path

# phenix.mtz.dump file
# phenix.mtz.dump --show_column_data file