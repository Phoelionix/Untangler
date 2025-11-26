set -u
RESNAME=$1
out_folder=$2

libcheck <<EOF
Y
_MON $RESNAME
_LIST L
_FILE_O $out_folder/$RESNAME
_END
EOF

rm $out_folder/${RESNAME}_${RESNAME}.pdb
rm $out_folder/${RESNAME}_${RESNAME}.cif
rm $out_folder/${RESNAME}_${RESNAME}.ps
rm $out_folder/${RESNAME}.odb