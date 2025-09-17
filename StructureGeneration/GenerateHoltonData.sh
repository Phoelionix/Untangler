set -u 

modelPDBPath=$1

cd $(dirname "$0")


if [ ! -f "$modelPDBPath" ]; then
    echo $modelPDBPath
    echo "filenotfound"
    exit 0
fi


mkdir -p HoltonScores
mkdir -p HoltonOutputs
cd HoltonOutputs

handle=$(basename "$modelPDBPath" .pdb) 
#rm -f ${handle}_log.txt

../../Measures/untangle_score.csh ../$modelPDBPath get_individual_residue_scores=True > "${handle}_log.txt"



score=$(tail -1 "${handle}_log.txt" | head -1 | awk '{print $1}')

#rm -f ../HoltonScores/*${handle}_log.txt

# Exclude water swapped files TODO make optional argument to this script
if [[ $handle != *_WaSw ]]; then 
    cp "${handle}_log.txt" "../HoltonScores/${score}-${handle}_log.txt"
fi

