#BSUB -J "28_8_100t_bh"
#BSUB -o "/home/maherm/bs_out/28_8_100t_bh.%J.o"
#BSUB -e "/home/maherm/bs_out/28_8_100t_bh.%J.e"
#BSUB -R "rusage[mem=3]"
#BSUB -n "10"
#BSUB -W 24:00
#BSUB -R "select[type==CentOS7]"
#BSUB -R span[hosts=1]

SCRIPTFP="/home/maherm/scgenome/scgenome/scripts/simulate_BHC.py"
VENVFP="/home/maherm/scgenome/venv/bin/activate"

echo "Sourcing ${VENVFP}"
source $VENVFP

python $SCRIPTFP
