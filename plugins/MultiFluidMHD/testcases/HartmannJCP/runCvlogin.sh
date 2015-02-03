#!/bin/bash -l
#SBATCH --ntasks=64 # 64 core per node
#SBATCH --exclusive  #no other jobs allowed on the node while this job is running
#SBATCH --time=48:00:00
#SBATCH -J STOAVT205 #job name

# #SBATCH --nodelist=node[025-028]
# #SBATCH -o /scratch/alaguna/MRStudy3_Quarter
# #SBATCH --nodes=1  # 64 cores per node

# settings
# NTASKS -> overall number of processors to be reserved are specified in --ntasks above
CASEDIR=$PWD
NCOREPERNODE=64
NPROCPERNODE=64

module load vki/coolfluid/openmpi/1.6.5  2> /dev/null

# change directory
cd $CASEDIR || exit 1
echo "ENTERING DIR: $PWD"

# prepare a hostfile and give some info
HOSTFILE=$CASEDIR/hosts.$SLURM_JOB_ID
rm -fv $HOSTFILE
# SORT!!!!!!!!!!!
srun hostname -s  | sort | awk "NR % ${NCOREPERNODE} > (( ${NCOREPERNODE} -1 )) - $NPROCPERNODE" > $HOSTFILE
NCORES=$( cat $HOSTFILE | wc -l )

# some system info
echo "NUMCORES: $NCORES"
echo "NODELIST: " $( srun hostname -s | sort | awk "NR % ${NCOREPERNODE} > (( ${NCOREPERNODE} -1 )) - 1" )
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
ulimit -a

# and run the stuff
# --mca mpi_paffinity_alone 1 \
echo "-----x-x-x- RUN -x-x-x-----"
export OPENBLAS_NUM_THREADS=1
export CASEFILE="./hartmannAirHalfTube.CFcase"

ln -s /home/alaguna/CF2/OPENMPI/optim/src/Solver/coolfluid-solver* ./ 
#EXE="$( mpirun -bind-to-core -report-bindings -np $NCORES -hostfile $HOSTFILE /home/alaguna/CF2/OPENMPI/optim/src/Solver/coolfluid-solver --scase $CASEFILE > output )"
EXE="$( mpirun -bind-to-core -report-bindings -np $NCORES -hostfile $HOSTFILE ./coolfluid-solver --scase $CASEFILE > output )"
echo $EXE
$EXE


