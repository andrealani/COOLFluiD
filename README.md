COOLFluiD
=========

The object-oriented HPC platform for CFD, plasma and multi-physics simulations whose development started in 2002 at the Von Karman Institute for Fluid Dynamics is finally open. 

#########################################################################
#                                                                       #
# This is some information to get started with COOLFluiD                #
# (more detailed info will be available on the website, once restored). #
#                                                                       #
#########################################################################

In order to install COOLFluiD you first need to install:

mpi      (e.g. openmpi, mpich, mvapich)  
petsc    (3.3 or 3.4.2)
boost    (1.42 or 1.54)
curl     (optional)
parmetis (4.0.3) 
 
Once you have done that, you must set the correct path inside the coolfluid.conf file (see examples in tools/conf/).

./prepare.pl --config-file=CONF_FILE --mods-update        (download all required modules)
./prepare.pl --config-file=CONF_FILE --build=BUILD_MODE   (configures build system and creates makefiles)

Note that if CONF_FILE=coolfluid.conf (default) you can skip "--config-file=CONF_FILE".

If you don't have errors, you should be ready to compile now:

cd basebuild_dir/BUILD_MODE   (see what you specified in CONF_FILE for "basebuild_dir")

make -jN

When the compilation is finished, your executable is:

basebuild_dir/BUILD_MODE/src/Solver/coolfluid_solver 

In order o run, you need to have also the 

basebuild_dir/BUILD_MODE/src/Solver/coolfluid_solver.xml

inside the same directory. It is convenient to just create soft links to 
coolfluid_solver* in the directory TESTDIR where collfluid needs to be run :

cd TESTDIR
ln -sf basebuild_dir/BUILD_MODE/src/Solver/coolfluid_solver* .

To run:

mpirun -np NPROC ./coolfluid_solver --scase TESTCASE_PATH/testcase.CFcase

TESTCASE_PATH must be "./" if the .CFcase file is in the current directory.
