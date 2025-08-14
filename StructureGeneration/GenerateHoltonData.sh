


modelPDBPath=$1

cd $(dirname "$0")

echo $(dirname "$0")

if [ ! -f $modelPDBPath ]; then
    echo $modelPDBPath
    echo "filenotfound"
    exit 0
fi


mkdir -p HoltonOutputs
cd HoltonOutputs

../../Measures/untangle_score.csh ../$modelPDBPath get_individual_residue_scores=True





# chmod a+x ../../UnTangle/untangle_score.csh
# ../..//UnTangle/untangle_score.csh manual_built.pdb debug=1 > manual_score.log
# tar czvf send_to_holton.tgz tempfile* manual*