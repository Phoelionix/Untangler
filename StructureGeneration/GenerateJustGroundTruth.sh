set -u

pdbID=$1
numConformations=3
shake_before_geo_optim=0.1
useCDL=False
#pdbID="6QIY" # 9MIZ not workin?

#TODO make directory data and data/tmp_refinement




# mapfile -t FIELDS < <(
#   awk -F', ' '{ for (i=1; i<=NF; i++) print $i }' "$(dirname "$0")/pdb_IDs.csv"
# )

#for pdbID in 2QEV 1S4R 4KUA 1S1Z 1S4S 4I38 2JLH 3JVI 4KOU 5E56 4ZMD 2G1K 8GEV 2AIF 3T8R 7RZB; do 
cd $(dirname "$0")/data

# mapfile -t FIELDS < <(
#   awk -F', ' '{ for (i=1; i<=NF; i++) print $i }' "$(dirname "$0")/pdb_IDs.csv"
# )
#for pdbID in "${FIELDS[@]}"; do
    bash ../FetchPDB.sh $pdbID
    bash ../GenerateConformationsGroundTruth.sh $pdbID $numConformations $shake_before_geo_optim $useCDL
    #bash ../GenerateConformationsForInitialRefine.sh $pdbID $numConformations 
    # bash ../ReadySet.sh $pdbID $numConformations 
    # python3.9 ../CombineStructuresToEnsemble.py ${pdbID}_ground_truth ${pdbID}_conf_1H.pdb ${pdbID}_conf_2H.pdb
    #python3.9 ../CombineStructuresToEnsemble.py ${pdbID}_initial_model ${pdbID}_starting_1H.pdb ${pdbID}_starting_2H.pdb
    
    structures=$(eval echo ${pdbID}_conf_{1..$numConformations}.pdb)

    python3.9 ../CombineStructuresToEnsemble.py ${pdbID}_simTruth-$numConformations $structures
    bash ../GenerateScatteringData.sh ${pdbID}_simTruth-$numConformations

    cd ../
    bash GenerateHoltonData.sh data/${pdbID}_simTruth-$numConformations.pdb
    # bash ../GenerateRefinement.sh $pdbID
    #python3.9 ../CreateStartingEnsemblePair.py $pdbID.pdb

    ####. ../GenerateInitialRefine.sh $pdbID # $numConformations
#done;