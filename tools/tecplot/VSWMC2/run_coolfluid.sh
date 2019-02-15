#!/bin/bash
#
#PBS -A lp_vswmc
# 
##PBS -l nodes=1:ppn=20:ivybridge
#PBS -l nodes=6:ppn=20:ivybridge
# 
### Regular queue
##PBS -l walltime=02:00:00
#PBS -l walltime=24:00:00
### Debug queue
##PBS -l qos=debugging
# 
#PBS -j oe
#
#PBS -N COOLFluiD_SolarWind
#
#PBS -m abe
#
#PBS -M andrey.kochanov@kuleuven.be

#module load bzip2/1.0.6-intel-2014a
source switch_to_2015a
module purge
module load cURL/7.43.0-foss-2015a
module load PETSc/3.6.3-foss-2015a-Python-2.7.9

cd $RUN_PATH
# example:
#RUN_PATH=/data/leuven/310/vsc31047/vswmc_runs/coolfluid/test1/
#cd /data/leuven/310/vsc31047/vswmc_runs/coolfluid/test1/

#CONFIG_FILE=${1:-storm.CFcase}
#LOG_FILE=${2:-SolarWind_log.txt}

ln -s /data/leuven/310/vsc31047/VSWMC_MODEL/COOLFluiD/OPTIM/optim/apps/Solver/coolfluid-solver* .

#mpirun -np 120 ./coolfluid-solver --scase ./storm.CFcase >& SolarWind_log.txt
#mpirun -np 20 ./coolfluid-solver --scase ./"$CONFIG_FILE" &> ./"$LOG_FILE"
mpirun -np 120 ./coolfluid-solver --scase ./"$CONFIG_FILE" &> ./"$LOG_FILE"

source switch_to_2014a
module purge
#module load Python/2.7.6-foss-2014a
module load K40c/2014a
module load Tecplot/2017R3

python CFplot_final_2.py $PWD/RESULTS storm  &>> ./"$LOG_FILE"

echo EOF >> ./"$LOG_FILE"

