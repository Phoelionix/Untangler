set -u

handle=$1



ground_truth=${handle}.pdb

phenix.fmodel $ground_truth high_resolution=1 type=real #add_sigmas=True type=real scale=0.03
mv $ground_truth.mtz ${handle}-sf.mtz
#phenix.reflection_file_converter ${pdbID}_conf_${i}-sf.mtz --sca=${pdbID}_ground_truth_scatter.sca THIS DOESNT WORK FOR SOME REASONNNN. Is it sigmas? Holton said ignores sigmas. Weird. 


#     awk 'NR>=10 && NR<=15' file1.txt > file2.txt


#--write_mtz_intensities

