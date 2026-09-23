#/bin/bash
#OAR --name femtic_annecy_inversion
#OAR --project mt-geothest
#OAR -l /nodes=1,walltime=40:0:00
#OAR -p total_mem=256
#OAR --notify mail:charroyj@univ-smb.fr

ulimit -s unlimited

source /applis/site/nix.sh
nix-env --switch-profile $NIX_USER_PROFILE_DIR/ced_intel
source $INTEL_ONEAPI/setvars.sh

#project=annecy

cd /home/charroyj/FemTic/Software/femticPy/projects/random_start/median/

mpirun -np 4 --machinefile $OAR_NODE_FILE -launcher-exec "oarsh" /home/charroyj/bin/femtic

#mpirun -np `cat $OAR_FILE_NODES|wc -l` --machinefile $OAR_NODE_FILE -launcher-exec "oarsh" /home/charroyj/bin/ifemtic.x

