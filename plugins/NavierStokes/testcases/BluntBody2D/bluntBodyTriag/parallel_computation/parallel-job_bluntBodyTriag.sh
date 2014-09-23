#!/bin/bash
# author Andrea Lani
# author Giuseppe Angelini
# author Thomas Wuilbaut
# author Tiago Quintino

uname -a
cd $HOME/workspace/COOLFluiD/optim/

# Get number of nodes, processors to $NPROC
export NPROC=$(wc -l < $PBS_NODEFILE)
export COOLFLUID_PROG="mpirun -np $NPROC ./coolfluid"
export COOLFLUID_DIR="$HOME/workspace/COOLFluiD/"
export COOLFLUID_BUILD_DIR="$HOME/workspace/COOLFluiD/optim/"

# Creation of machines_cyl file
MACHINE_FILE = machines_.$(date +%H%M%S%d%m%Y)
cat $PBS_NODEFILE > $MACHINE_FILE
echo agclust1.local >> $MACHINE_FILE

lamwipe $MACHINE_FILE
lamboot $MACHINE_FILE

date
./run-coolfluid.sh ../../../testcases/TestGambit/first_evaluation_meeting/Blunt_Body_2D/bluntBodyTriag/bluntBodyLTE.CFcase
date

lamhalt
lamwipe $MACHINE_FILE
