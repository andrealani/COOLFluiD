#!/bin/bash -l
#SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_vswmc"
#SBATCH --cluster=wice
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
#SBATCH --nodes="3"
#SBATCH --ntasks-per-node="90"
#SBATCH --ntasks="270"
#SBATCH --time="4:00:00"

PMIX_MCA_gds=hash
module load CMake/3.20.1-GCCcore-10.3.0
module load Boost/1.76.0-GCC-10.3.0
module load ParMETIS/4.0.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a

mpirun ./coolfluid-solver --scase /scratch/leuven/357/vsc35794/Coronal_modelling/TimeDependent/COCONUT_Map/map_steady_fullMHD.CFcase 


