#!/bin/bash

source switch_to_2015a
module purge
module load cURL/7.43.0-foss-2015a
module load PETSc/3.6.3-foss-2015a-Python-2.7.9

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
cd "$SCRIPTPATH"

CONFIG_FILE=${1:-storm_low.CFcase}
LOG_FILE=${2:-SolarWind_log.txt}

ln -s /data/leuven/310/vsc31047/VSWMC_MODEL/COOLFluiD/OPTIM/optim/apps/Solver/coolfluid-solver* .

#mpirun -np 20 ./coolfluid-solver --scase ./storm_low.CFcase >& SolarWind_log.txt
mpirun -np 20 ./coolfluid-solver --scase ./"$CONFIG_FILE" &> ./"$LOG_FILE"

source switch_to_2014a
module purge
#module load Python/2.7.6-foss-2014a
module load K40c/2014a
module load Tecplot/2017R3

python CFplot_final_2.py $PWD/RESULTS storm &>> ./"$LOG_FILE"

echo EOF >> ./"$LOG_FILE"

