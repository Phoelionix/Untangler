pdbID=$1

numConformations=$2

# Readyset - Add H, etc.
# TODO this should be a separate step. 
for i in $(seq 1 $numConformations); # todo change to query for all files matching pattern of: ${pdbID}_conf_*.cif
do
    coordFile=${pdbID}_conf_${i}H.pdb
    refineFile=${pdbID}_starting_${i}H.pdb

    phenix.ready_set ${pdbID}_conf_$i.pdb
    phenix.ready_set ${pdbID}_starting_$i.pdb
    mv ${pdbID}_conf_${i}.updated.pdb $coordFile
    mv ${pdbID}_starting_${i}.updated.pdb $refineFile
    rm ${pdbID}_conf_${i}.updated.cif
    rm ${pdbID}_starting_${i}.updated.cif

done;

# this removes the anisotropy.
echo "WARNING assuming 2 configs"
tag=$pdbID
python3.9 ../CombineStructuresToEnsemble.py $tag ${pdbID}_conf_1H.pdb ${pdbID}_conf_2H.pdb 
sleep 10
echo ${tag}_ground_truth.pdb
