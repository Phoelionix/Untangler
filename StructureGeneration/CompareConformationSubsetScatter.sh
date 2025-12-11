set -u


pdb_path1=$1
pdb_path2=$2
altlocs1=$3 # Which altlocs should be considered in scattering?

# Option to specify altlocs for pdb_path2
if [ -z "${4+xxx}" ]; then
 altlocs2=$altlocs1
else
 altlocs2=$4
fi

# Generate R factor for conformation structure factors


pdb_file1=${pdb_path1##*/}
pdb_handle1=${pdb_file1%.*}
pdb_file2=${pdb_path2##*/}
pdb_handle2=${pdb_file2%.*}

subset1=$(dirname "$0")/output/$pdb_handle1-$altlocs1.pdb
subset2=$(dirname "$0")/output/$pdb_handle2-$altlocs2.pdb

rm -f $subset1 $subset2 
bash $(dirname "$0")/make_altloc_subset.sh $pdb_path1 $altlocs1 $subset1 > /dev/null
bash $(dirname "$0")/make_altloc_subset.sh $pdb_path2 $altlocs2 $subset2 > /dev/null

mtz_subset2=$(dirname "$0")/output/$pdb_handle2-$altlocs2.mtz
rm -f $mtz_subset2
bash $(dirname "$0")/GenerateScatteringData.sh $subset2 $mtz_subset2 

#bash $(dirname "$0")/GetRData.sh $subset1 $mtz_subset2
bash $(dirname "$0")/PrintRData.sh $subset1 $mtz_subset2