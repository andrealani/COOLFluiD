#########################################################################################
# Generic OPTIONS
#########################################################################################

# user option to add assertions
OPTION ( CF_ENABLE_ASSERTIONS "Enable code assertions" ON )
IF ( NOT CF_ENABLE_ASSERTIONS)
  ADD_DEFINITIONS(-DNDEBUG)
ENDIF ()

# precision real numbers
IF ( NOT CF_PRECISION_SINGLE )
OPTION(CF_PRECISION_SINGLE       "Real numbers have single precision"       OFF )
ENDIF()
IF ( NOT CF_PRECISION_LONG_DOUBLE )
OPTION(CF_PRECISION_LONG_DOUBLE   "Real numbers have long double precision" OFF )
ENDIF()
IF ( NOT CF_PRECISION_SINGLE AND NOT CF_PRECISION_LONG_DOUBLE )
OPTION(CF_PRECISION_DOUBLE   "Real numbers have double precision"           ON )
ENDIF()

# user option to add tracing
OPTION(CF_ENABLE_TRACE 	"Enable tracing code"  ON)
IF( NOT CF_ENABLE_TRACE)
  ADD_DEFINITIONS(-DCF_NO_TRACE)
ENDIF()

# user option to add logging
OPTION(CF_ENABLE_LOGALL 	"Enable logging via CFLog facility" ON)
IF( NOT CF_ENABLE_LOGALL)
  ADD_DEFINITIONS(-DCF_NO_LOG)
ENDIF()

# user option to add debug logging
OPTION(CF_ENABLE_LOGDEBUG 	"Enable debug logging via CFLog facility" ON)
IF( NOT CF_ENABLE_LOGDEBUG)
  ADD_DEFINITIONS(-DCF_NO_DEBUG_LOG)
ENDIF()

# user option to enable debug macros
OPTION(CF_ENABLE_DEBUG_MACROS 	"Enable debug macros"                 ON)
IF( NOT CF_ENABLE_DEBUG_MACROS)
  ADD_DEFINITIONS(-DCF_NO_DEBUG_MACROS)
ENDIF()

# user option to search other dirs for plugins
SET ( CF_EXTRA_SEARCH_DIRS "" CACHE STRING "Full paths to extra dirs to be searched for plugin modules which maybe out of source." )

OPTION ( CF_ENABLE_MPI                "Enable MPI compilation"                  ON   )
OPTION ( CF_ENABLE_DOCS               "Enable build of documentation"           ON   )
OPTION ( CF_ENABLE_EXPLICIT_TEMPLATES "Enable explicit template instantiation"  ON   )
OPTION ( CF_ENABLE_GROWARRAY          "Enable GrowArray usage"                  ON   )
OPTION ( CF_ENABLE_INTERNAL_DEPS      "Enable internal dependencies between libraries"  ON   )
OPTION ( CF_ENABLE_AUTOMATIC_UPDATE_MODULES  "Enable automatic subversion update of the plugins" OFF  )
OPTION ( CF_ENABLE_TESTCASES          "Enable checking testcases from CMake system" ON )
OPTION ( CF_ENABLE_UNITTESTS          "Enable creation of unit tests"            OFF )
OPTION ( CF_ENABLE_WARNINGS           "Enable lots of warnings while compiling"  ON )
OPTION ( CF_ENABLE_STDASSERT          "Enable standard assert() functions "  ON )

OPTION ( CF_ENABLE_PARALLEL_VERBOSE   "Enable extra output in the parallel interface" OFF  )
OPTION ( CF_ENABLE_PARALLEL_DEBUG     "Enable debug code on the parallel interface"  OFF  )

OPTION ( CF_CMAKE_LIST_PLUGINS             "CMake lists the plugins"                 OFF  )

SET    ( CF_TESTCASES_NCPUS "2" CACHE STRING "Number of CPUs to ue in parallel run" )

# if user disables MPI we overwrite the CF_HAVE_MPI variable
IF ( CF_ENABLE_MPI )
  IF ( CF_MPI_AVAILABLE )
    SET ( CF_HAVE_MPI 1 CACHE BOOL "User enabled MPI [FOUND]" )
  ELSE()
    SET ( CF_HAVE_MPI 0 CACHE BOOL "User enabled MPI [NOT-FOUND]" )
  ENDIF()
ELSE()
  SET ( CF_HAVE_MPI 0 CACHE BOOL "User disabled MPI" )
ENDIF ()

#########################################################################################
# PROFILING OPTIONS
#########################################################################################

# user option to add system depedent profiling
OPTION ( CF_ENABLE_PROFILING    "Enable code profiling"                 OFF )

IF(CF_ENABLE_PROFILING)

  ###########################
  # GNU gprof
  IF(CF_PROFILER_TOOL MATCHES gprof)
    IF(UNIX AND CMAKE_COMPILER_IS_GNUCC)
      SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -pg" )
      SET(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -pg" )
    ELSE(UNIX AND CMAKE_COMPILER_IS_GNUCC)
      LOG("User selected profiler [gprof] must be used with GCC compiler")
      SET( CF_PROFILER_TOOL     NOTFOUND )
    ENDIF()
  ENDIF()

  ###########################
  # google-perftools
  IF(CF_PROFILER_TOOL MATCHES google-perftools)

    FIND_PACKAGE(GooglePerftools)

    IF(CF_HAVE_GOOGLE_PERFTOOLS)
      LINK_LIBRARIES(${GOOGLE_PERFTOOLS_LIBRARIES})
    ELSE(CF_HAVE_GOOGLE_PERFTOOLS)
      LOG("User selected profiler [google-pertools] could not be found")
      SET( CF_PROFILER_TOOL     NOTFOUND )
    ENDIF()

  ENDIF()

ENDIF()

#########################################################################################
# STATIC BUILD OPTIONS
#########################################################################################

# user option to static build
OPTION ( CF_ENABLE_STATIC       "Enable static building"                OFF)

IF ( CF_ENABLE_STATIC )

  SET(BUILD_SHARED_LIBS OFF)

    IF(UNIX)
      # LINUX
      IF("${CMAKE_SYSTEM}" MATCHES Linux)
        SET(CMAKE_CXX_LINK_EXECUTABLE
        "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> -Wl,-whole-archive <LINK_LIBRARIES> -Wl,-no-whole-archive")
      ENDIF("${CMAKE_SYSTEM}" MATCHES Linux)
      # SGI IRIX
      IF("${CMAKE_SYSTEM}" MATCHES IRIX)
        SET(CMAKE_CXX_LINK_EXECUTABLE
        "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> -Wl,-all <LINK_LIBRARIES> -Wl,-notall")
      ENDIF("${CMAKE_SYSTEM}" MATCHES IRIX)
      # On Darwin:
      #  -all_load $convenience
      IF("${CMAKE_SYSTEM}" MATCHES Darwin)
        SET(CMAKE_CXX_LINK_EXECUTABLE
        "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> -all_load <LINK_LIBRARIES>")
      ENDIF("${CMAKE_SYSTEM}" MATCHES Darwin)
      # On Solaris 2:
      #   -z allextract $convenience -z defaultextract
    ENDIF(UNIX)

ELSE()

  SET( BUILD_SHARED_LIBS ON )

ENDIF()




