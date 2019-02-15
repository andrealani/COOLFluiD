#!/bin/bash

#module load tecplot/360_2009
source switch_to_2014a
module purge
module load K40c/2014a
module load Tecplot/2017R3

#cd /software/alternate/coolfluid/VSWMC2_CASES/TestSteady
#cd /data/leuven/310/vsc31047/vswmc_runs/coolfluid/test1/

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
cd "$SCRIPTPATH"

preplot RESULTS/Steady06042000Storm_Final.plt RESULTS/Steady06042000Storm_Final.plt.bin
#preplot RESULTS/Steady06042000Storm_Final.surf.plt RESULTS/Steady06042000Storm_Final.surf.plt.bin 
#export pathtoscript=/data/leuven/310/vsc31047/vswmc_runs/coolfluid/test1

OLD="/software/alternate/coolfluid/VSWMC2_CASES/TestSteady"

sed -i "s|"$OLD"|"$SCRIPTPATH"|g" CFplot_final.mcr

tec360 -b -p CFplot_final.mcr

# the following images are generated: 
#RESULTS/Bz_iter_0.png  (Bz) 
#RESULTS/p_iter_0.png   (pressure) 
#RESULTS/rho_iter_0.png (density)

#gimp *.png&
