# Add some deviations. Minimize geometry energy.

pdbID=$1

numConformations=$2

minimizeGeometry=$3

tag=$4

shake=$5


for i in $(seq 1 $numConformations);
do
   phenix.pdbtools $pdbID.pdb shake=$shake
        
    

    if [ "$minimizeGeometry" = true ] ; then
        phenix.geometry_minimization ${pdbID}_modified.pdb
        mv ${pdbID}_modified_minimized.pdb ${pdbID}_${tag}_$i.pdb
        rm ${pdbID}_modified_minimized.geo
        rm ${pdbID}_modified_minimized.cif
    else

        #python3.9 ../cif2pdb.py ${pdbID}_modified.cif ${pdbID}_modified.pdb $referenceCellGeomPdbPath
        mv ${pdbID}_modified.pdb ${pdbID}_${tag}_$i.pdb
        #rm ${pdbID}_modified.cif
    fi
    rm ${pdbID}_modified.pdb


    

   
     

done;