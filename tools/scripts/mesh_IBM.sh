# @ job_name = Example
# @ comment = "Example Job"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ wall_clock_limit = 00:30:00
# @ notification = error
# @ notify_user = lani@vki.ac.be
# @ job_type = bluegene
# @ bg_size = 60
# @ queue

rm -fr core.*
#module load petsc/3.7.6_basic_O3simd_C++_int8 >& /dev/null
#module load parmetis/4.0.3 >& /dev/null
#module load boost/1.47.0 >& /dev/null

#export LD_LIBRARY_PATH=/homeb/prpc19/prpc1900/local/petsc_3.6.3/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH is $LD_LIBRARY_PATH"

#export OMP_NUM_THREADS=64
#runjob --ranks-per-node 32 --exp-env OMP_NUM_THREADS --exe ./example_exe

export CFEXEC=/homeb/prpc19/prpc1900/trunk_LANI/LANI_opt2/bgoptim2/apps/Solver/coolfluid-solver
echo "Running the executable: $CFEXEC"

# --log 650 corresponds to CFLog(VERBOSE, ...) in the code
runjob --np 960 --ranks-per-node 16 : $CFEXEC --scase ./extrude.CFcase
# --log 650
# --testEnv true
