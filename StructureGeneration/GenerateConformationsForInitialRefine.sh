pdbID=$1
numConformations=$2
useCDL=$3
minimizeGeometry=true
shake=0.1

bash ../GenerateConformations.sh $pdbID $numConformations $minimizeGeometry starting $shake $useCDL