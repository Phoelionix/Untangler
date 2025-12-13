set -u 

modelPDBPath=$1; shift 1

cd $(dirname "$0")

cdl='true'
keep_hydrogens='true'
Rwork_override='false'
Rfree_override='false'
while getopts ":F:W:CH" flag; do
 case $flag in
    F) Rfree=$OPTARG
       Rfree_override='true'
    ;;
    W) Rwork=$OPTARG
       Rwork_override='true'
    ;;
    C) cdl='false'
    ;;
    H) keep_hydrogens='false'
    ;;
   \?)
   echo INVALID FLAG
   exit 0
   ;;
 esac
done

wxray=1 #50 #1 #10

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

out_handle=$handle
if ! $keep_hydrogens; then 
    out_handle=${out_handle}_ignoreH
    args="${args} keep_hydrogens=$keep_hydrogens"
fi
if ! $cdl; then 
    out_handle=${out_handle}_noCDL
fi

args="${args} outprefix=$out_handle"


if $Rwork_override; then
    args="${args} Rwork=$Rwork"
fi 
if $Rfree_override; then
    args="${args} Rfree=$Rfree"
fi 


args="${args} wxray=$wxray cdl=$cdl"


#../../Measures/untangle_score.csh ../$modelPDBPath $args > "${out_handle}_log.txt" 
../../Measures/untangle_score_weighted.csh ../$modelPDBPath $args > "${out_handle}_log.txt" 



score=$(tail -1 "${out_handle}_log.txt" | head -1 | awk '{print $1}')

#rm -f ../HoltonScores/*${out_handle}_log.txt

# Exclude water swapped files TODO make optional argument to this script
if [[ $handle != *_WaSw ]]; then 
    cp "${out_handle}_log.txt" "../HoltonScores/${score}-${out_handle}_log.txt"
fi

