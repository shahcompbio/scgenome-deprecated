# Script for submitting multiple jobs doing BHC, each job on a different
# dataset
#BSUB -J tyler_sets_bhc[1-11]
#BSUB -o "/home/maherm/bs_out/tyler_sets_bhc.%J.o"
#BSUB -e "/home/maherm/bs_out/tyler_sets_bhc.%J.e"
#BSUB -R "rusage[mem=8]"
#BSUB -n "1"
#BSUB -W 48:00
#BSUB -R "select[type==CentOS7]"
#BSUB -R span[hosts=1]


SCRIPTFP="/home/maherm/scgenome/scgenome/scripts/do_bhc.py"
#SCRIPTFP="/home/maherm/test_arg.py"
VENVFP="/home/maherm/scgenome/venv/bin/activate"
ARGFP="/work/shah/maherm/tantalus/ado_bhc_subfile"
ARGLINE=`head -n $LSB_JOBINDEX $ARGFP | tail -n 1`

echo "Sourcing ${VENVFP}"
source $VENVFP

python $SCRIPTFP $ARGLINE
