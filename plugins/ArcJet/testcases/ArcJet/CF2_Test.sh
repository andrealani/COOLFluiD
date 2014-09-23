#PBS -S /bin/csh
#PBS -N CF2_Test
# This example uses the Westmere nodes
# User job can access ~22.5 GB of memory per Westmere node.
# A memory intensive job that needs more than ~1.9 GB
# per process should use less than 16 cores per node
# to allow more memory per MPI process. This example
# asks for 1 node and 2 MPI processes per node.
# This request implies 1x2 = 2 MPI processes for the job.
#PBS -l select=1:ncpus=2:mpiprocs=2:model=wes
#PBS -l walltime=00:05:00
#PBS -q devel 
#PBS -j oe
#PBS -W group_list=s1340
#PBS -m e

# Loading modules 

module load comp-intel/2013.5.192
module load mpi-sgi/mpt.2.08r7

# By default, PBS executes your job from your home directory.
# However, you can use the environment variable
# PBS_O_WORKDIR to change to the directory where
# you submitted your job.

cd $PBS_O_WORKDIR

# use of dplace to pin processes to processors may improve performance

mpiexec -np 2 ./coolfluid-solver --scase ./arcjet_phi_Mut1_LTE_Test.CFcase > output
# Example of use dplace: mpiexec dplace -s1 -c 4-11 ./grinder < run_input > output

# It is a good practice to write stderr and stdout to a file (ex: output)
# Otherwise, they will be written to the PBS stderr and stdout in /PBS/spool,
# which has limited amount  of space. When /PBS/spool is filled up, any job
# that tries to write to /PBS/spool will die.

# -end of script-
