set -u
RESNAME=$1
out_folder=$2

# pdb_file=${pdb_path##*/}
# pdb_handle=${pdb_file%.*}

# libcheck <<EOF
# Y
# _FILE_PDB $pdb_path
# _MON *
# _LIST L
# _END
# EOF
#_FILE_O $out_folder/$pdb_handle


# _MON  j1
# _FILEL  j1min.lib
# _fileo j1new


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