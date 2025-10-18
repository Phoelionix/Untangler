set -u 

modelPDBPath=$1

cd $(dirname "$0")


wxray=1 #10

#auto_ignore_geometries='true'
auto_ignore_geometries='false'


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

args=""
if $auto_ignore_geometries; then
    args="${args} overridefile=${handle}_potential_overrides.txt"
fi

args="${args} wxray=$wxray"


#../../Measures/untangle_score.csh ../$modelPDBPath $args > "${handle}_log.txt" 
../../Measures/untangle_score_weighted.csh ../$modelPDBPath $args > "${handle}_log.txt" 



score=$(tail -1 "${handle}_log.txt" | head -1 | awk '{print $1}')

#rm -f ../HoltonScores/*${handle}_log.txt

# Exclude water swapped files TODO make optional argument to this script
if [[ $handle != *_WaSw ]]; then 
    cp "${handle}_log.txt" "../HoltonScores/${score}-${handle}_log.txt"
fi

