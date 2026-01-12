
sigf_mtz=~/Untangler/data/refme.mtz
target_mtz=~/Untangler/data/1aho_simTruth-4-sf.mtz
out=~/Untangler/StructureGeneration/sigf_add_test.mtz

FP=FMODEL


Rflag=FreeR_flag
SIGFP=SIGF


if [ -f $out ]; then 
mv $out $out#
fi

sftools << EOF 
read $target_mtz
delete  col $Rflag
read $sigf_mtz col $SIGFP $Rflag
absent col $FP if col $SIGFP absent
absent col $Rflag if col $FP absent
absent col $SIGFP if col $FP absent
select col $SIGFP = PRESENT
purge nodata yes
select all
write $out col $FP $SIGFP $Rflag
quit
y
EOF


