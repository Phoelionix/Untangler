# bash Refine.sh ../data/4et8.pdb refined -d ../data/lys_galli_aug_real.mtz -n 5 -s 3 -gzt

# TODO folder with file containing commands used, shell file, and .eff file. 

set -u

xyz_path=$1 
hkl_path=$2; shift 2 

xyz_file=${xyz_path##*/}
xyz_handle=${xyz_file%.*}

trials=0

dampA=0.02
dampB=0.05

####

out_handle_override='false'
unrestrained='false'
calc_wE='false'
refine_water_occupancy='false'

while getopts ":o:n:A:B:uwW" flag; do
 case $flag in
    o) out_handle=$OPTARG
       out_handle_override='true'
    ;;
    n) trials=$OPTARG
    ;;
    A) dampA=$OPTARG
    ;;
    B) dampB=$OPTARG
    ;;
    u) unrestrained='true'
    ;;
    w) calc_wE='true'
    ;;
    W) refine_water_occupancy='true'
    ;;
   \?)
   echo INVALID FLAG
   exit 0
   ;;
 esac
done


if ! $out_handle_override; then
  out_handle=${xyz_handle}_refined
fi 


cd $(dirname "$0")

working_folder=$out_handle

mkdir -p tmp_refinement


mkdir -p tmp_refinement/$working_folder

cp $xyz_path tmp_refinement/$working_folder/initial_model.pdb
cp $hkl_path tmp_refinement/$working_folder/xray.mtz


cd tmp_refinement/$working_folder


logs_path="../../../output/refine_logs"
mkdir -p $logs_path

######Options######

cat > refmac_opts.txt  << EOF 
make hydr Y
make hout Y
EOF

# Appending...

cat >> refmac_opts.txt << EOF
damp $dampA $dampA $dampB
EOF

if $unrestrained; then
cat << EOF >> refmac_opts.txt
weight matrix 99
REFI TYPE UNREstrained
EOF
else
cat >> refmac_opts.txt << EOF 
weight matrix 0.5
EOF
fi

# cat > refmac_opts.txt  << EOF 
# damp 0.02 0.02 0.05
# make hydr Y
# make hout Y
# weight matrix 0.5
# EOF



if $refine_water_occupancy; then
setup_occupancy_waters.com initial_model.pdb refmac_opts.txt
fi 
##############################

#########
refmacout=refmacout_minRfree.pdb

if [ -f $refmacout ]; then
  mv $refmacout ${refmacout}#
fi

if [ "$trials" -eq "0" ]; then 
  converge_refmac.com "initial_model.pdb xray.mtz" > $logs_path/${out_handle}.log
else
  converge_refmac.com "initial_model.pdb xray.mtz trials=$trials" > $logs_path/${out_handle}.log
fi

if [ ! -f $refmacout ]; then
  echo Could not find output file $refmacout
  exit 0
fi

cp $refmacout ../../../output/${out_handle}.pdb  

cd ../.. 

if $calc_wE; then
  cd ../StructureGeneration
  bash GenerateHoltonData.sh ../output/${out_handle}.pdb  > HoltonScores/${out_handle}.log
fi









######## setup_occupancy_waters.com #######
# ##
# #! /bin/tcsh -f
# #
# #
# #

# if(! -e "$2") then
#     echo "usage: $0 refmacin.pdb refmac_opts.txt"
#     exit 9
# endif

# grep HOH $1 >! occme.pdb
# refmac_occupancy_setup.com occme.pdb >> $2