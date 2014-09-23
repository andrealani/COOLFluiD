#######################################################################
# find MPI compiler or libraries

# try to compile an mpi program to check if compiler is already mpi
CHECK_CXX_SOURCE_COMPILES (
"#include <mpi.h>
 #include <iostream>
 int main(int argc, char* argv[])
 {
   MPI_Init(&argc, &argv); int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Finalize();
   return 0;
 }"
CF_MPI_COMPILER_AVAILABLE)

IF ( CF_MPI_COMPILER_AVAILABLE )

  LOG ( "Already using MPI C++ compiler, no need to search for MPI libraries" )
  SET ( CF_MPI_AVAILABLE 1 CACHE BOOL "Found MPI library" )

ELSE  ()

  LOG ( "No MPI C++ compiler was set. Searching for MPI libraries..." )
  FIND_PACKAGE(MPI)
  IF ( NOT CF_MPI_AVAILABLE)
    LOG ( "No MPI libraries found." )
    LOG ( "WARNING: user asked for MPI but MPI is deactivated" )
  ENDIF( NOT CF_MPI_AVAILABLE)

ENDIF ()

#######################################################################
# find mpirun

FIND_PROGRAM(CF_MPIRUN_PROGRAM mpirun
            PATHS ${MPI_HOME}/bin $ENV{MPI_HOME}/bin
            PATH_SUFFIXES mpi/bin
            DOC "mpirun program"
            NO_DEFAULT_PATH )
FIND_PROGRAM(CF_MPIRUN_PROGRAM mpirun
            PATH_SUFFIXES mpi/bin
            DOC "mpirun program" )

IF    (NOT CF_MPIRUN_PROGRAM)
  SET (CF_MPIRUN_PROGRAM mpirun CACHE STRING "mpirun program set by default")
ENDIF ()
MARK_AS_ADVANCED ( CF_MPIRUN_PROGRAM )
LOGFILE ( "  CF_MPIRUN_PROGRAM : [${CF_MPIRUN_PROGRAM}]" )
