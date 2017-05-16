##############################################################################
# COOLFluiD cmake build system
##############################################################################
# TODO
#  * packaging with the build system
#  * check precompiled headers
##############################################################################
# Command line options for cmake:
# -DCMAKE_BUILD_TYPE=DEBUG     will set the build type
# -DCF_SKIP_FORTRAN=:BOOL=ON   will not test for fortran language
# -DCMAKE_SKIP_RPATH:BOOL=ON   will remove rpath
# -G"Visual Studio 9 2008"     will generate for WIN32
# -DCMAKE_INSTALL_PREFIX=DIR   will specify where to install
##############################################################################

##############################################################################
# cmake behavior
##############################################################################

CMAKE_MINIMUM_REQUIRED( VERSION 2.8.3 FATAL_ERROR )

# Tell the CMake makefile generator to not have rules depend on themselves.
# This causes extra rebuilds when the include path changes from turning a library on or off
SET(CMAKE_SKIP_RULE_DEPENDENCY 1)
# allow relaxed loop and if constructions
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# check for cmake policies introduced in cmake 2.6
if(COMMAND cmake_policy)

  # check 'cmake --help-policy CMP000*'

  if ( CMAKE_MINOR_VERSION GREATER 5 ) # cmake > 2.6

    # cmake >= 2.6.0

    # policy ignores CMAKE_BACKWARDS_COMPATIBILITY keyword
    cmake_policy(SET CMP0001 NEW)
    # policy checks unique global names for targets
    cmake_policy(SET CMP0002 NEW)
    # policy prefers libraries with full path when linking libraries
    cmake_policy(SET CMP0003 NEW)
    # policy gives error on leading and trailing spaces in library names
    cmake_policy(SET CMP0004 NEW)
    # policy tries to escape the -D compile definitions
    cmake_policy(SET CMP0005 NEW)
    # policy asks for bundle destination in MACOSX
    cmake_policy(SET CMP0006 NEW)
    # policy ignores empty elements on a list
    cmake_policy(SET CMP0007 OLD)

    # cmake >= 2.6.1
    if (CMAKE_PATCH_VERSION GREATER 0)
      # policy passes full path library names to the generator
      cmake_policy(SET CMP0008 NEW)
    endif()

    # cmake > 2.6.2
    if (CMAKE_PATCH_VERSION GREATER 1)
      # policy glob_recurse does not follow symbolic links
      cmake_policy(SET CMP0009 NEW)
    endif()

    # cmake >= 2.6.3
    if (CMAKE_PATCH_VERSION GREATER 2)
      # bad variable references treated as errors
      cmake_policy(SET CMP0010 NEW)
      # included scripts do cmake_policy push and pop
      cmake_policy(SET CMP0011 NEW)
    endif()
  endif()

  if ( CMAKE_MINOR_VERSION GREATER 7 ) # cmake > 2.8

    # cmake >= 2.8.0

    # if() recognizes numbers and boolean constants
    cmake_policy(SET CMP0012 NEW)
    # Duplicate binary directories are not allowed
    cmake_policy(SET CMP0013 NEW)
    # Input directories must have CMakeLists.txt
    cmake_policy(SET CMP0014 NEW)

  endif()

endif(COMMAND cmake_policy)

##############################################################################
# project definition and versioning
##############################################################################

PROJECT (COOLFluiD C CXX )

# include personal options
INCLUDE(${COOLFluiD_BINARY_DIR}/coolfluid.cmake OPTIONAL)

# option to skip fortran language
IF( NOT CF_SKIP_FORTRAN )
  ENABLE_LANGUAGE(Fortran)
ENDIF()

SET ( CF_SOURCE_EXTENSIONS cxx cpp cc c f for f77 f90 )

# option to add CUDA support
IF ( CF_ENABLE_CUDA )
  # the following is needed by .../share/cmake-2.8/Modules/CMakeDetermineCUDACompiler.cmake 
  SET ( CMAKE_PLATFORM_ROOT_BIN "${CMAKE_BINARY_DIR}/CMakeFiles/${CMAKE_VERSION}")
  MESSAGE ("CMAKE_PLATFORM_ROOT_BIN : ${CMAKE_PLATFORM_ROOT_BIN}")
  ENABLE_LANGUAGE(CUDA)
  LIST ( APPEND CF_SOURCE_EXTENSIONS cu )
ENDIF()

SET ( CF_HEADER_EXTENSIONS ci hh h hpp )
SET ( CF_FILE_EXTENSIONS ${CF_SOURCE_EXTENSIONS} ${CF_HEADER_EXTENSIONS} )

SET (COOLFLUID_VERSION_MAJOR    2013)
SET (COOLFLUID_VERSION_MINOR    9)
SET (COOLFLUID_VERSION "${COOLFLUID_VERSION_MAJOR}.${COOLFLUID_VERSION_MINOR}" )

SET (CF_KERNEL_VERSION_MAJOR    2)
SET (CF_KERNEL_VERSION_MINOR    5)
SET (CF_KERNEL_VERSION_MICRO    0)
SET (CF_KERNEL_VERSION          "${CF_KERNEL_VERSION_MAJOR}.${CF_KERNEL_VERSION_MINOR}.${CF_KERNEL_VERSION_MICRO}" )

##############################################################################
# cmake setup
##############################################################################

# disallow in-source build
IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "COOLFluiD requires an out of source build.\nPlease create a separate build directory and run 'cmake path_to_coolfluid [options]' there.")
ENDIF()

# include other configurations
SET(CF_CMAKE_DIR "${COOLFluiD_SOURCE_DIR}/cmake")
SET(CMAKE_MODULE_PATH "${CF_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

INCLUDE( DefineMacros )   # define macros
INCLUDE( CTest )          # include testing support

# recreate the logfile
SET  ( PROJECT_LOG_FILE ${PROJECT_BINARY_DIR}/CMakeLogInfo.txt )
FILE ( WRITE ${PROJECT_LOG_FILE} "COOLFluiD Cmake Log file\n")

# Append the library version information to the library target properties.
# A parent project may set its own properties and/or may block this.
IF(NOT CF_NO_LIBRARY_VERSION)
  SET(COOLFluiD_LIBRARY_PROPERTIES ${COOLFluiD_LIBRARY_PROPERTIES}
    VERSION "${COOLFluiD_VERSION}"
    SOVERSION "${COOLFluiD_VERSION_MAJOR}.${COOLFluiD_VERSION_MINOR}"
    )
ENDIF()

MESSAGE ( STATUS "CMAKE_INSTALL_PREFIX [${CMAKE_INSTALL_PREFIX}]"  ) 

IF ( NOT DEFINED CF_INSTALL_BIN_DIR )
  SET(CF_INSTALL_BIN_DIR bin)
ENDIF ()

IF ( NOT DEFINED CF_INSTALL_LIB_DIR )
  SET(CF_INSTALL_LIB_DIR lib)
ENDIF ()

IF ( NOT DEFINED CF_INSTALL_INCLUDE_DIR )
  SET(CF_INSTALL_INCLUDE_DIR include/coolfluid)
ENDIF ()

IF ( NOT DEFINED CF_INSTALL_DATA_DIR )
  SET(CF_INSTALL_DATA_DIR data)
ENDIF ()

# use RPATHs, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# create the dso directory for shared libraries
SET ( COOLFluiD_DSO_DIR ${COOLFluiD_BINARY_DIR}/dso )
FILE ( MAKE_DIRECTORY ${COOLFluiD_DSO_DIR} )

# name the generated config files
SET ( CF_CONFIG_HH_FILENAME coolfluid_config.h CACHE INTERNAL "name of configuration file" )
SET ( CF_SVN_HH_FILENAME    coolfluid_svnversion.hh CACHE INTERNAL "name of svn version file" )

##############################################################################
# configuration options
##############################################################################

# set the default build to be RelWithDebInfo
IF ( NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF()

INCLUDE ( CheckOperatingSystem  )   # check for operating system
INCLUDE ( CheckCompilerFeatures )   # check compiler features

INCLUDE ( DefineGlobalOptions   )   # add user global options
INCLUDE ( DefineBuildModes      )   # create extra build modes
INCLUDE ( DefineBuildRules      )   # defined default compilation flags and linking rules per architecture

##############################################################################
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
LOG ( "General Configuration ")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
##############################################################################

ADD_SUBDIRECTORY ( cmake )
ADD_SUBDIRECTORY ( doc )

##############################################################################

IF (CF_ENABLE_STATIC)
  SET(CF_HAVE_ALLSTATIC 1 CACHE BOOL "All libraries compiled with static linking")
ELSE()
  SET(CF_HAVE_ALLSTATIC 0 CACHE BOOL "All libraries compiled with dynamic linking")
ENDIF()

IF ( CF_ENABLE_IBMSTATIC )
   SET (CF_HAVE_IBMSTATIC 1 CACHE BOOL "Compiling on IBM machine with static linking")                                         
ELSE()  
   SET (CF_HAVE_IBMSTATIC 0 CACHE BOOL "Compiling on IBM machine with dynamic linking")
ENDIF()

IF ( CF_ENABLE_SINGLEEXEC ) 
   SET (CF_HAVE_SINGLE_EXEC 1 CACHE BOOL "Compile a single executable, no libraries")
ELSE()
   SET (CF_HAVE_SINGLE_EXEC 0 CACHE BOOL "Don't compile a single executable, no libraries")
ENDIF()

# finding boost (essential)

# find boost package (essential)
# SET ( Boost_DEBUG 1 )  # to debug boost search
SET ( Boost_USE_STATIC_LIBS ${CF_ENABLE_STATIC} )
SET ( Boost_USE_MULTITHREAD ON  )
# find based on minimal version defined below
SET ( Boost_FIND_VERSION        ON   )
SET ( Boost_FIND_VERSION_MAJOR  "1"  )
SET ( Boost_FIND_VERSION_MINOR  "42" )
SET ( Boost_FIND_VERSION_PATCH  "0"  )
# older cmakes dont have these versions
SET ( Boost_ADDITIONAL_VERSIONS  "1.60" "1.60.0" "1.59" "1.59.0" "1.55" "1.55.0" "1.54" "1.54.0" "1.47" "1.47.0" "1.42" "1.42.0" "1.41" "1.41.0"  "1.40" "1.40.0" "1.39" "1.39.0" )
# components to search for

#IF ( CF_ENABLE_IBMSHARED OR CF_ENABLE_IBMSTATIC ) 
#LIST ( APPEND CF_Boost_COMPONENTS filesystem system regex unit_test_framework )
#ELSE()
LIST ( APPEND CF_Boost_COMPONENTS thread filesystem system regex unit_test_framework )
#ENDIF()
FIND_PACKAGE( Boost COMPONENTS ${CF_Boost_COMPONENTS} )

LOG ( "Boost include path [${Boost_INCLUDE_DIR}]" )
LOG ( "Boost lib version  [${Boost_LIB_VERSION}]" )
LOG ( "Boost libraries    [${Boost_LIBRARIES}]" )

IF (${Boost_LIB_VERSION} STREQUAL "1_60" )
ADD_DEFINITIONS  ( -DCF_HAVE_BOOST_1_60 )
ELSEIF (${Boost_LIB_VERSION} STREQUAL "1_59" )
ADD_DEFINITIONS  ( -DCF_HAVE_BOOST_1_59 )
ELSEIF (${Boost_LIB_VERSION} STREQUAL "1_55" )
ADD_DEFINITIONS  ( -DCF_HAVE_BOOST_1_55 )
ELSEIF (${Boost_LIB_VERSION} STREQUAL "1_54" )
ADD_DEFINITIONS  ( -DCF_HAVE_BOOST_1_54 )
ELSEIF (${Boost_LIB_VERSION} STREQUAL "1_47")
ADD_DEFINITIONS  ( -DCF_HAVE_BOOST_1_47 )
ELSEIF (${Boost_LIB_VERSION} STREQUAL "1_42")
ADD_DEFINITIONS  ( -DCF_HAVE_BOOST_1_42 )
ENDIF()

# if not found give more information
IF ( NOT Boost_FOUND )
  LOG ( ${Boost_ERROR_REASON} )
  MESSAGE ( FATAL_ERROR "Boost is required to compile COOLFluiD Kernel" )
ENDIF()
# add boost include path
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIR} )
SET ( CF_Boost_LIBRARIES ${Boost_LIBRARIES} )

IF ( CF_ENABLE_SINGLEEXEC )
# remove debug version of Boost libraries if present
# we assume that the release version is available
 LIST ( REMOVE_ITEM CF_Boost_LIBRARIES debug )
 LIST ( REMOVE_ITEM CF_Boost_LIBRARIES optimized )
 FOREACH (ABoostLIB ${Boost_LIBRARIES})
  IF (${ABoostLIB} MATCHES "-d-")
   LIST ( REMOVE_ITEM CF_Boost_LIBRARIES ${ABoostLIB} )
  ENDIF()
 ENDFOREACH(ABoostLIB)
ENDIF()

##############################################################################
# find non essential packages

# our find macros
#FIND_PACKAGE(BlasLapack)      # search for Blas Lapack support
FIND_PACKAGE(Metis)           # serial domain decomposition
FIND_PACKAGE(Parmetis)        # parallel domain decomposition
#ADD_DEFINITIONS  ( -DIDXTYPEWIDTH=64 )   # 64 bit
#ADD_DEFINITIONS  ( -DREALTYPEWIDTH=64 )  # 64 bit

IF ( CF_ENABLE_CURL )
  FIND_PACKAGE(Curl)          # curl downloads files on the fly
ENDIF()
FIND_PACKAGE(Valgrind)        # valgrind for profiling and memmory leak detection
FIND_PACKAGE(GooglePerftools) # Google PerfTools for profiling

IF ( CF_ENABLE_GSL )
  FIND_PACKAGE(GSL)             # GNU Scientific Library
ENDIF() 

# MPI 
IF (CF_ENABLE_MPI)  
FIND_PACKAGE(MPI)
# if mpi was found add it to the include path if needed
IF   (CF_HAVE_MPI AND NOT CF_HAVE_MPI_COMPILER)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_DIR})
  ADD_DEFINITIONS  ( -DCF_HAVE_MPI )
ENDIF ()
ENDIF ()

IF (CF_ENABLE_CUDA AND CF_ENABLE_VIENNACL)
 ADD_DEFINITIONS  ( -DCF_HAVE_VIENNACL )
ENDIF()

# mutation++
IF ( CF_ENABLE_MUTATIONPP )
FIND_PACKAGE(Mutationpp)    # mutation++ library
# define CF_HAVE_MUTATIONPP only if enabled to do so
ADD_DEFINITIONS  ( -DCF_HAVE_MUTATIONPP )
ENDIF()

# PLATO
IF ( CF_ENABLE_PLATO )
FIND_PACKAGE(PLATO)    # PLATO library
# define CF_HAVE_PLATO only if enabled to do so
ADD_DEFINITIONS  ( -DCF_HAVE_PLATO )
ENDIF()

# CUDA
IF ( CF_ENABLE_CUDA )
FIND_PACKAGE(CUDA)  # CUDA environment
INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIR})
ENDIF()

# cmake find macros

FIND_PACKAGE(ZLIB)          # file compression support
LOG ( "ZLIB_FOUND: [${ZLIB_FOUND}]" )
IF ( ZLIB_FOUND )
	LOG ( "  ZLIB_INCLUDE_DIRS: [${ZLIB_INCLUDE_DIRS}]" )
	LOG ( "  ZLIB_LIBRARIES:    [${ZLIB_LIBRARIES}]" )
ENDIF()


# AL: this should not be here !!!!!!
IF ( NOT CF_ENABLE_SINGLEEXEC )
FIND_PACKAGE (BZip2)
LOG ( "BZIP2_FOUND: [${BZIP2_FOUND}]" )
IF ( BZIP2_FOUND )
	LOG ( "  BZIP2_INCLUDE_DIR:  [${BZIP2_INCLUDE_DIR}]" )
	LOG ( "  BZIP2_LIBRARIES:    [${BZIP2_LIBRARIES}]" )
	LOG ( "  BZIP2_DEFINITIONS:  [${BZIP2_DEFINITIONS}]" )
	LOG ( "  BZIP2_NEED_PREFIX:  [${BZIP2_NEED_PREFIX}]" )
ENDIF()
ENDIF()


#include the generated config file
ADD_DEFINITIONS  ( -DCF_HAVE_CONFIG_H )
INCLUDE_DIRECTORIES ( ${COOLFluiD_BINARY_DIR} )

# AL: SVN is not checked anymore for portability reasons
# check for subversion support
# INCLUDE(CheckSvnVersion)

# enable all kind of tests
IF ( CF_ENABLE_TESTCASES )
  SET (CF_ENABLE_UNIT_TESTS ON) 
  SET (CF_ENABLE_UNIT_CASES ON)    
  SET (CF_ENABLE_PERFORMANCE_TESTS ON)     

  SET (CF_ENABLE_PERFORMANCE_CASES ON)
ENDIF()

##############################################################################
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
LOG ( "Kernel Configuration ")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
##############################################################################

# continue to kernel modules
INCLUDE_DIRECTORIES( ${COOLFluiD_SOURCE_DIR}/src )
# include in compilation path the src dir
ADD_SUBDIRECTORY ( src )

#MESSAGE("##### coolfluid-solver_files = ${coolfluid-solver_files}")

##############################################################################
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
LOG ( "Plugin's Configuration ")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
##############################################################################

# check existence of Mutation* libraries and set corresponding environmental variables
IF (EXISTS "${COOLFluiD_SOURCE_DIR}/plugins/Mutation") 
  SET(CF_HAVE_MUTATION1 1 CACHE BOOL "Found Mutation library")
  LOG ( "${COOLFluiD_SOURCE_DIR}/plugins/Mutation exists" )
ELSE()
  SET(CF_HAVE_MUTATION1 0 CACHE BOOL "Not found Mutation library")
ENDIF()

IF (EXISTS "${COOLFluiD_SOURCE_DIR}/plugins/Mutation2.0") 
  SET(CF_HAVE_MUTATION2 1 CACHE BOOL "Found Mutation2 library")
  LOG ( "${COOLFluiD_SOURCE_DIR}/plugins/Mutation2.0 exists" )
ELSE()
  SET(CF_HAVE_MUTATION2 0 CACHE BOOL "Not found Mutation2 library")
ENDIF()

IF (EXISTS "${COOLFluiD_SOURCE_DIR}/plugins/Mutation2.0.0") 
  SET(CF_HAVE_MUTATION2OLD 1 CACHE BOOL "Found Mutation2OLD library")
  LOG ( "${COOLFluiD_SOURCE_DIR}/plugins/Mutation2.0.0 exists" )
ELSE()
  SET(CF_HAVE_MUTATION2OLD 0 CACHE BOOL "Not found Mutation2OLD library")
ENDIF()

# local plugins
INCLUDE_DIRECTORIES( ${COOLFluiD_SOURCE_DIR}/plugins )
ADD_SUBDIRECTORY ( plugins )

# find extra plugins
FOREACH( EXDIR ${CF_EXTRA_SEARCH_DIRS} )
  # include base extra dir in compilation search path
  INCLUDE_DIRECTORIES(${EXDIR})

  # include plugin configuration
  FILE(GLOB_RECURSE PLUGIN_CMAKE_CFGS "${EXDIR}/*/*.cmake")
  LIST(SORT PLUGIN_CMAKE_CFGS)
  FOREACH( ACFG ${PLUGIN_CMAKE_CFGS} )
    INCLUDE(${ACFG})
  ENDFOREACH(  ACFG ${PLUGIN_CMAKE_CFGS} )

  LOG ( "" )
  LOG ( "---------------------------------------------------------" )
  LOG ( "EXTRA MODULE DIR [${EXDIR}]")
  LOG ( "---------------------------------------------------------" )

  # find modules in each search dir
  FILE(GLOB_RECURSE EXTRA_PLUGIN_MODULES "${EXDIR}/*/CMakeLists.txt")
  LIST(SORT EXTRA_PLUGIN_MODULES)

  FOREACH( ADIR ${EXTRA_PLUGIN_MODULES} )
    STRING (REGEX REPLACE "(/)*CMakeLists.txt" "" MODDIR ${ADIR})
    IF (IS_DIRECTORY ${MODDIR} )
      # include each module dir in compilation search path
      INCLUDE_DIRECTORIES(${MODDIR})
      FILE(RELATIVE_PATH FINALDIR ${COOLFluiD_SOURCE_DIR} ${MODDIR})
      LIST ( APPEND CF_MODULES_LIST ${FINALDIR} )
      LIST ( APPEND CF_EXTRA_MODULES_LIST ${FINALDIR})
    ENDIF (IS_DIRECTORY  ${MODDIR} )
  ENDFOREACH( ADIR ${EXTRA_PLUGIN_MODULES} )
ENDFOREACH( EXDIR ${CF_EXTRA_SEARCH_DIRS} )

# include modules directories
FOREACH( amodule ${CF_EXTRA_MODULES_LIST} )
    LOG ( "\#---------------------------------------------------------" )
    LOG ( "\# MODULE [${amodule}]")
    ADD_SUBDIRECTORY ( ${amodule} extra )
ENDFOREACH( amodule ${CF_EXTRA_MODULES_LIST} )
LOG ( "---------------------------------------------------------" )

# applications directory (Solver & co)
ADD_SUBDIRECTORY ( apps )

#config file generation (at the end, so everything is detected)
ADD_SUBDIRECTORY ( config )

# install the headers
INSTALL(FILES ${COOLFluiD_BINARY_DIR}/${CF_CONFIG_HH_FILENAME}
              ${COOLFluiD_BINARY_DIR}/${CF_SVN_HH_FILENAME}
       DESTINATION
              ${CF_INSTALL_INCLUDE_DIR})

##############################################################################
# print summary
##############################################################################

LOG ( "")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
LOG ( "+ COOLFluiD configuration summary ")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
LOG ( " COOLFluiD Release     : [${COOLFLUID_VERSION}]")
LOG ( " COOLFluiD Kernel      : [${CF_KERNEL_VERSION}]")
LOG ( " COOLFluiD SVN version : [${coolfluid_svnversion}]")
LOG ( "---------------------------------------------------------" )
LOG ( " CMake Version         : [${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}]")
LOG ( " CMake Generator       : [${CMAKE_GENERATOR}]")
LOG ( " Build Type            : [${CMAKE_BUILD_TYPE}]")
LOG ( " Operating System      : [${CMAKE_SYSTEM}]")
LOG ( " Memory address        : [${CF_OS_BITS} bits]")
LOG ( "---------------------------------------------------------" )
LOG ( " C compiler            : [${CMAKE_C_COMPILER}]")
LOG ( " C flags               : [${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}]")
LOG ( " C++ compiler          : [${CMAKE_CXX_COMPILER}]")
LOG ( " C++ flags             : [${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}]")
LOG ( " C++ link flags        : [${CMAKE_CXX_LINK_FLAGS}]")
IF( NOT CF_SKIP_FORTRAN )
LOG ( " Fortran compiler      : [${CMAKE_Fortran_COMPILER}]")
LOG ( " Fortran flags         : [${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}]")
ENDIF()
IF ( CF_ENABLE_CUDA )
LOG ( " CUDA compiler      : [${CMAKE_CUDA_COMPILER}]")
LOG ( " CUDA flags         : [${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE}}]")
LOG ( " CUDA malloc        : [${CF_HAVE_CUDA_MALLOC}]")
ENDIF()
LOG ( " Shared Linker flags   : [${CMAKE_SHARED_LINKER_FLAGS}]")
LOG ( "---------------------------------------------------------" )
LOG ( " MPI                   : [${CF_HAVE_MPI}]" )
IF(CF_HAVE_MPI)
LOG ( "    CF_MPIRUN_PROGRAM  : [${CF_MPIRUN_PROGRAM}]" )
LOG ( "    MPI_INCLUDE_DIR    : [${MPI_INCLUDE_DIR}]")
LOG ( "    MPI_LIBRARIES      : [${MPI_LIBRARIES}]")
ENDIF()
LOG ( "---------------------------------------------------------" )
LOG ( " Check internal deps   : [${CF_ENABLE_INTERNAL_DEPS}]")
LOG ( " Log all               : [${CF_ENABLE_LOGALL}]")
LOG ( " Log debug             : [${CF_ENABLE_LOGDEBUG}]")
LOG ( " Assertions            : [${CF_ENABLE_ASSERTIONS}]")
LOG ( " Tracing               : [${CF_ENABLE_TRACE}]")
LOG ( " Static libs           : [${CF_ENABLE_STATIC}]")
LOG ( " Profiling             : [${CF_ENABLE_PROFILING}]")
LOG ( " long int              : [${CF_HAVE_LONG}]")
LOG ( " long long int         : [${CF_HAVE_LLONG}]")
LOG ( " CURL enabled          : [${CF_ENABLE_CURL}]")
LOG ( " CUDA enabled          : [${CF_ENABLE_CUDA}]")
LOG ( " BOOST libs            : [${CF_Boost_LIBRARIES}]") 
IF(CF_ENABLE_PROFILING)
LOG ( "    Profiler           : [${CF_PROFILER_TOOL}]")
ENDIF()
LOG ( " Explicit Templates    : [${CF_HAVE_CXX_EXPLICIT_TEMPLATES}]")
LOG ( "---------------------------------------------------------" )
LOG ( " Kernel libs           : [${CF_KERNEL_LIBS}]")
LOG ( " Static plugin libs    : [${CF_KERNEL_STATIC_LIBS}]")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
LOG ( "")

# if file with orphan files exists remove it
SET ( ORPHAN_FILE "${COOLFluiD_BINARY_DIR}/OrphanFiles.txt" )
IF ( EXISTS ${ORPHAN_FILE} ) 
	FILE ( REMOVE ${ORPHAN_FILE} )
ENDIF()
# if orphan files were found, put the list on the file
LIST(LENGTH CF_ORPHAN_FILES CF_LENGTH_ORPHAN_FILES)
IF (CF_LENGTH_ORPHAN_FILES)
  LOG ( " !!! WARNING !!! Orphan files were found during the configuration.")
  LOG ( " !!! WARNING !!! Check full list in file ${ORPHAN_FILE} ")
  FOREACH( AFILE ${CF_ORPHAN_FILES} )
    FILE ( APPEND ${ORPHAN_FILE} "${AFILE}\n" )
  ENDFOREACH()
ENDIF()

DUMP_VARIABLES()        # dump relevant variables in the log file
INCLUDE(PrepareCPack)   # packaging instructions

