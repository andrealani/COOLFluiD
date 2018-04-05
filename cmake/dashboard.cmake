##########################################################################
# General Settings
################################################################

# directories
set(CF_ROOT /data1/banyai/codes/coolfluid2/nightly)
set(CF_DEPS /data1/banyai/codes/coolfluid2/local)

# number of cores to use
set(CF_NBPROCS 4)

# user and machine names 
set(CTEST_SITE "$ENV{USER}")
set(CTEST_BUILD_NAME "$ENV{HOSTNAME}")

##########################################################################
# Advanced Settings
#
# Requires that you know what you doing
##########################################################################

# dashboard model 
set(MODEL Nightly)

# subdirs to put in
set(CTEST_SOURCE_DIRECTORY "${CF_ROOT}/source")
set(CTEST_BINARY_DIRECTORY "${CF_ROOT}/build")

# ctest config
set(CF_BUILD_TYPE "ReleaseWithDebInfo")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# miscellaneous helpers
set(CTEST_BUILD_FLAGS "-j${CF_NBPROCS}")
set(CTEST_UPDATE_COMMAND "svn")
set(CTEST_CMAKE_COMMAND "${CF_DEPS}/bin/cmake")

# miscellaneous locations/files
set(CF_SOURCE_SVN "https://coolfluidsrv.vki.ac.be/svn/coolfluid/Sources/Kernel/trunk")
set(CF_CONF_FILE "coolfluid.conf.ctest")

# configure options
set(CF_CONFIG_OPTIONS
  -DCMAKE_C_COMPILER=${CF_DEPS}/bin/mpicc 
  -DCMAKE_CXX_COMPILER=${CF_DEPS}/bin/mpicxx 
  -DBOOST_ROOT=${CF_DEPS} 
  -DCMAKE_INSTALL_PREFIX=${CF_ROOT}/install 
  -DCMAKE_BUILD_TYPE=${CF_BUILD_TYPE} 
  -DCF_ENABLE_DOCS=ON 
  -DCF_ENABLE_UNIT_TESTS=ON 
  -DCF_ENABLE_PERFORMANCE_TESTS=ON 
  -DCF_ENABLE_UNIT_CASES=ON 
  -DCF_ENABLE_PERFORMANCE_CASES=ON 
  -DCF_TESTING_NB_PROCS=${CF_NBPROCS} 
  -DPETSC_DIR=${CF_DEPS}/petsc 
)
message(${CF_CONFIG_OPTIONS})

# switches
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY ON)
set(CF_SUBMIT_TO_DASHBOARD_SERVER ON)
set(CF_PERFORM_UPDATE   ON)
set(CF_PERFORM_CONFIG   ON)
set(CF_PERFORM_BUILD    ON)
set(CF_PERFORM_UTESTS   ON)
set(CF_PERFORM_PTESTS   OFF)
set(CF_PERFORM_UCASES   ON)
set(CF_PERFORM_PCASES   OFF)
set(CF_PERFORM_COVERAGE OFF)

##########################################################################
# Pre 
##########################################################################
message("START")

# start
if ( CTEST_START_WITH_EMPTY_BINARY_DIRECTORY )
  ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
endif()

# enrich path
set(ENV{PATH} "${CF_DEPS}/bin:$ENV{PATH}")
set(ENV{LD_LIBRARY_PATH} "${CF_DEPS}/lib:$ENV{LD_LIBRARY_PATH}")

# flags for coverage
if (CF_PERFORM_COVERAGE)
  SET(CTEST_COVERAGE_COMMAND "gcov")
  list(APPEND CF_CONFIG_OPTIONS "-DCMAKE_C_FLAGS=--coverage" )
  list(APPEND CF_CONFIG_OPTIONS "-DCMAKE_CXX_FLAGS=--coverage" )
  list(APPEND CF_CONFIG_OPTIONS "-DCMAKE_C_LINK_FLAGS=--coverage" )
  list(APPEND CF_CONFIG_OPTIONS "-DCMAKE_CXX_LINK_FLAGS=--coverage" )
endif()

# very careful check with nbprocs because make -j without number is a guaranteed machine crash
if (${CF_NBPROCS} STREQUAL "" )
  message(FATAL_ERROR "CF_NBPROCS is not an  integer.")
endif()
string( REGEX REPLACE "[0-9]" "" test_nbprocs "${CF_NBPROCS}")
if ( NOT ( "${test_nbprocs}" STREQUAL "" ) )
  message(FATAL_ERROR "CF_NBPROCS is not an  integer.")
endif()

##########################################################################
# Checkout & initiate
##########################################################################

# since plugin checkout is hidden, checkout is going to be reported into framework
# et this point ctestconfig.txt doesn't exist, so using execute_process to perform initial checkout
if ( NOT EXISTS ${CTEST_SOURCE_DIRECTORY} )
message("CLEAN CHECKOUT")

  set(MODEL Experimental)
  configure_file(${CF_CONF_FILE} ${CTEST_SOURCE_DIRECTORY}/coolfluid.conf @ONLY)
  execute_process(COMMAND ${CTEST_UPDATE_COMMAND} co ${CF_SOURCE_SVN} ${CTEST_SOURCE_DIRECTORY} )
  execute_process(COMMAND ${CTEST_SOURCE_DIRECTORY}/prepare.pl --mods-update WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})

endif()

#initiate
set(START_TIME ${CTEST_ELAPSED_TIME})
include("${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake")
set_property(GLOBAL PROPERTY Label Everything)
ctest_start (${MODEL})

##########################################################################
# Update
##########################################################################

# collect plugin directories into a list
set(CF_PLUGIN_DIRS "")
file( GLOB CF_PLUGIN_DIRS "${CTEST_SOURCE_DIRECTORY}/plugins/[a-z,A-Z,0-9]*" )
foreach(idir ${CF_PLUGIN_DIRS})
  if( NOT (IS_DIRECTORY ${idir}) )
    list(REMOVE_ITEM CF_PLUGIN_DIRS ${idir})
  endif()
endforeach()

# update is all-in-one 
# this could be done per-plugin, but ctest_update doesn't like it
if( CF_PERFORM_UPDATE )
message("UPDATE")

  set(updatedirs "${CTEST_SOURCE_DIRECTORY}")
  foreach(iplug ${CF_PLUGIN_DIRS})
    set(updatedirs "${updatedirs} ${iplug}")
  endforeach()
  set(CTEST_CHECKOUT_COMMAND "${CTEST_UPDATE_COMMAND} up ${updatedirs}")
  ctest_update(APPEND)
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Update)
  endif()

endif()


##########################################################################
# Config
##########################################################################

# configure is all-in-one 
if( CF_PERFORM_CONFIG )
message("CONFIG")

  # perform config
  set(confopts "")
  foreach(iopt ${CF_CONFIG_OPTIONS})
    set(confopts "${confopts} ${iopt}")
  endforeach()

  set(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND} ${confopts} ${CTEST_SOURCE_DIRECTORY}")
  ctest_configure(APPEND)
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Configure)
  endif()

  # write project.xml
  file(WRITE    "${CTEST_BINARY_DIRECTORY}/Project.xml" "<Project name=\"coolfluid2\">\n")
  file(APPEND   "${CTEST_BINARY_DIRECTORY}/Project.xml" "  <SubProject name=\"Framework\">\n")
  file(APPEND   "${CTEST_BINARY_DIRECTORY}/Project.xml" "  </SubProject>\n")
  foreach(iplug ${CF_PLUGIN_DIRS})
    set(iplugname ${iplug})
    string(REGEX REPLACE "${CTEST_SOURCE_DIRECTORY}/plugins/" "" iplugname ${iplugname})
    file(APPEND "${CTEST_BINARY_DIRECTORY}/Project.xml" "  <SubProject name=\"${iplugname}\">\n")
    #file(APPEND "${CTEST_BINARY_DIRECTORY}/Project.xml" "    <Dependency name=\"Framework\">\n")
    #file(APPEND "${CTEST_BINARY_DIRECTORY}/Project.xml" "    </Dependency>\n")
    file(APPEND "${CTEST_BINARY_DIRECTORY}/Project.xml" "  </SubProject>\n")
  endforeach()
  file(APPEND   "${CTEST_BINARY_DIRECTORY}/Project.xml" "</Project>\n")

endif()
if (CF_SUBMIT_TO_DASHBOARD_SERVER)
  ctest_submit(FILES "${CTEST_BINARY_DIRECTORY}/Project.xml")
endif()

##########################################################################
# Build
##########################################################################

# compilation 
# this could be done per-plugin, if dependencies were known
if( CF_PERFORM_BUILD )
message("BUILD")

  ctest_build(BUILD ${CTEST_BINARY_DIRECTORY} APPEND)
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Build)
  endif()

endif()

##########################################################################
# Run tests
##########################################################################

if( CF_PERFORM_UTESTS OR CF_PERFORM_PTESTS OR CF_PERFORM_UCASES OR CF_PERFORM_PCASES )
  message("TEST")
endif()

# Framework's tests
set_property(GLOBAL PROPERTY SubProject Framework)
set_property(GLOBAL PROPERTY Label Framework)

if( CF_PERFORM_UTESTS )
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/src" APPEND INCLUDE "test-unit")
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Test)
  endif()
endif()
if( CF_PERFORM_PTESTS )
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/app" APPEND INCLUDE "test-perf")
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Test)
  endif()
endif()
if( CF_PERFORM_UCASES )
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/src" APPEND INCLUDE "case-unit")
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Test)
  endif()
endif()
if( CF_PERFORM_PCASES )
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/apps" APPEND INCLUDE "case-perf")
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Test)
  endif()
endif

# Testing plugins
foreach(iplug ${CF_PLUGIN_DIRS})

  set(iplugname ${iplug})
  string(REGEX REPLACE "${CTEST_SOURCE_DIRECTORY}/plugins/" "" iplugname ${iplugname})
  set_property(GLOBAL PROPERTY SubProject ${iplugname})
  set_property(GLOBAL PROPERTY Label ${iplugname})

  if( CF_PERFORM_UTESTS )
    ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/plugins/${iplugname}" APPEND INCLUDE "test-unit")
    if (CF_SUBMIT_TO_DASHBOARD_SERVER)
      ctest_submit(PARTS Test)
    endif()
  endif()
  if( CF_PERFORM_PTESTS )
    ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/plugins/${iplugname}" APPEND INCLUDE "test-perf")
    if (CF_SUBMIT_TO_DASHBOARD_SERVER)
      ctest_submit(PARTS Test)
    endif()
  endif()
  if( CF_PERFORM_UCASES )
    ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/plugins/${iplugname}" APPEND INCLUDE "case-unit")
    if (CF_SUBMIT_TO_DASHBOARD_SERVER)
      ctest_submit(PARTS Test)
    endif()
  endif()
  if( CF_PERFORM_PCASES )
    ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/plugins/${iplugname}" APPEND INCLUDE "case-perf")
    if (CF_SUBMIT_TO_DASHBOARD_SERVER)
      ctest_submit(PARTS Test)
    endif()
  endif()

endforeach()

##########################################################################
# Coverage
##########################################################################

if( CF_PERFORM_COVERAGE )
message("COVERAGE")

  # coverage of framework
  set_property(GLOBAL PROPERTY SubProject Framework)
  set_property(GLOBAL PROPERTY Label Framework)

  ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}/apps" APPEND)
 
  if (CF_SUBMIT_TO_DASHBOARD_SERVER)
    ctest_submit(PARTS Coverage)
  endif()

  # coverage of plugins
  foreach(iplug ${CF_PLUGIN_DIRS})

    set(iplugname ${iplug})
    string(REGEX REPLACE "${CTEST_SOURCE_DIRECTORY}/plugins/" "" iplugname ${iplugname})
    set_property(GLOBAL PROPERTY SubProject ${iplugname})
    set_property(GLOBAL PROPERTY Label ${iplugname})

    ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}/plugins/${iplugname}" APPEND)
 
    if (CF_SUBMIT_TO_DASHBOARD_SERVER)
      ctest_submit(PARTS Coverage)
    endif()

  endforeach()

endif()

##########################################################################
# Post 
##########################################################################
message("FINISH")
