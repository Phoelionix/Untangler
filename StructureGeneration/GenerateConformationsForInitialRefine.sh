pdbID=$1
numConformations=$2
minimizeGeometry=true
shake=0.1

bash ../GenerateConformations.sh $pdbID $numConformations $minimizeGeometry starting $shake