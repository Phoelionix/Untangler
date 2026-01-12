set -u

pdb_path=$1

pdb_file=${pdb_path##*/}
handle=${pdb_file%.*}

echo "YO"

# Option to specify out_path
if [ -z "${2+xxx}" ]; then
 out_path=${handle}-sf.mtz
else
 out_path=$2
fi


phenix.fmodel $pdb_path high_resolution=1 type=real #add_sigmas=True type=real scale=0.03
mv $pdb_file.mtz $out_path
#phenix.reflection_file_converter ${pdbID}_conf_${i}-sf.mtz --sca=${pdbID}_ground_truth_scatter.sca THIS DOESNT WORK FOR SOME REASONNNN. Is it sigmas? Holton said ignores sigmas. Weird. 


#     awk 'NR>=10 && NR<=15' file1.txt > file2.txt


#--write_mtz_intensities

