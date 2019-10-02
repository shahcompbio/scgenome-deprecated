#BSUB -J "ado_bhc"
#BSUB -o "/home/maherm/bs_out/ado_bhc.%J.o"
#BSUB -e "/home/maherm/bs_out/ado_bhc.%J.e"
#BSUB -R "rusage[mem=50]"
#BSUB -n "1"
#BSUB -W 24:00
#BSUB -R "select[type==CentOS7]"
#BSUB -R span[hosts=1]

SCRIPTFP="/home/maherm/scgenome/scgenome/scripts/do_bhc.py"
VENVFP="/home/maherm/scgenome/venv/bin/activate"

echo "Sourcing ${VENVFP}"
source $VENVFP

python $SCRIPTFP
