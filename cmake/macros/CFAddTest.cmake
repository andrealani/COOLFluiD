# Function to add a test.
#
# Mandatory keywords:
# - UTEST/PTEST
#      choose according to the test profile: unit-test, performance-test
#      give as value the name of the test
# - CPP
#      lists the sources and headers
#
# Optional keywords:
# - ARGUMENTS
#      a list of arguments to pass to the the test execution (default empty)
# - LIBS
#      a list of required libraries for this test to build (default empty)
# - MPI (default: mpirun is not called , unless global CMake constant makes all utests run with mpirun -np 1)
#      number of processors to use for call with mpirun (it can be a list of numbers)
#      the keyword "default" will take the global number CF_TESTING_NB_PROCS
# - DEPENDS
#      list of targets this test depends on (LIBS are automatically a dependency already)
#
# The two master switches turns the each type on and off
#   - CF_ENABLE_UNIT_TESTS
#   - CF_ENABLE_PERFORMANCE_TESTS
#
# After calling this function, the test is added to one of the following lists:
#   - CF_ENABLED_UTESTS
#   - CF_DISABLED_UTESTS
#   - CF_ENABLED_PTESTS
#   - CF_DISABLED_PTESTS
# The following variables will be set:
#   - ${TEST_TARGETNAME}_builds
#   - ${TEST_TARGETNAME}_dir
#   - ${TEST_TARGETNAME}_libs

function( cf_add_test )

#MESSAGE(INFO,"***TEST*************************************************************")
#MESSAGE(INFO,"***TEST*** CF_ENABLE_UNIT_TESTS= ${CF_ENABLE_UNIT_TESTS}")
#MESSAGE(INFO,"***TEST*** CF_ENABLE_PERFORMANCE_TESTS= ${CF_ENABLE_PERFORMANCE_TESTS}")

  set( single_value_args UTEST PTEST)
  set( multi_value_args  CPP ARGUMENTS MPI LIBS DEPENDS)
  set( _TEST_DIR ${CMAKE_CURRENT_BINARY_DIR} )

  # parse and complain if stg wrong with the arguments
  cmake_parse_arguments(_PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN})
  foreach(unparsed_arg ${_PAR_UNPARSED_ARGUMENTS})
    if( NOT ( (${unparsed_arg} STREQUAL "AND") OR (${unparsed_arg} STREQUAL "OR") OR (${unparsed_arg} STREQUAL "NOT") ) )
      list( APPEND _UNPARSED_ARGUMENTS ${unparsed_arg} )
    endif()
  endforeach()
  if(_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unknown keywords given to cf_add_test(): \"${_UNPARSED_ARGUMENTS}\"")
  endif()
  if( (NOT _PAR_UTEST) AND (NOT _PAR_PTEST))
    message(FATAL_ERROR "The call to cf_add_test() doesn't set the required \"UTEST/ATEST/PTEST test-name\" argument.")
  endif()

  # select unit/acceptance/performance tests
  if(_PAR_UTEST)
    set(_TEST_NAME ${_PAR_UTEST})
    set(_TEST_TARGETNAME "test-unit-${_PAR_UTEST}")
    set(_TEST_ENABLED ${CF_ENABLE_UNIT_TESTS})
  elseif(_PAR_PTEST)
    set(_TEST_NAME ${_PAR_PTEST})
    set(_TEST_TARGETNAME "test-perf-${_PAR_PTEST}")
    set(_TEST_ENABLED ${CF_ENABLE_PERFORMANCE_TESTS})
  endif()

#MESSAGE(INFO,"***TEST*** _TEST_NAME= ${_TEST_NAME}")
#MESSAGE(INFO,"***TEST*** _TEST_ENABLED= ${_TEST_ENABLED}")

  # option to build it or not (option is advanced and does not appear in the cmake gui)
  option( CF_BUILD_${_TEST_NAME} "Build the ${_TEST_PROFILE} [${_TEST_NAME}]" ${_TEST_ENABLED} )
  mark_as_advanced(CF_BUILD_${_TEST_NAME})

  # check if test will build
  set(_TEST_BUILDS ${CF_BUILD_${_TEST_NAME}})
  if( _PAR_MPI AND (NOT CF_HAVE_MPI) )
    set(_TEST_BUILDS FALSE)
  endif()
  if( (DEFINED _PAR_MPI) AND (NOT CF_HAVE_MPI) )
    set(_TEST_BUILDS FALSE)
  endif()

  # select mpi stuff
  set(_RUN_MPI FALSE)
  if (DEFINED _PAR_MPI)
    set(_RUN_MPI TRUE)
    if (NOT DEFINED CF_TESTING_NB_PROCS)
      set(CF_TESTING_NB_PROCS 2)
    endif()
  endif()

#MESSAGE(INFO,"***TEST*** _TEST_BUILDS= ${_TEST_BUILDS}")
#MESSAGE(INFO,"***TEST*** Boost_UNIT_TEST_FRAMEWORK_LIBRARY= ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
#MESSAGE(INFO,"***TEST*** Boost_UNIT_TEST_FRAMEWORK_FOUND= ${Boost_UNIT_TEST_FRAMEWORK_FOUND}")
#MESSAGE(INFO,"***TEST*** CF_HAVE_MPI= ${CF_HAVE_MPI}")
#MESSAGE(INFO,"***TEST*** MPI_LIBRARIES= ${MPI_LIBRARIES}")
#MESSAGE(INFO,"***TEST*** BOOSTUNITTEST_LIBRARY= ${BOOSTUNITTEST_LIBRARY}")

  # prepare
  set(_TEST_COMMAND ${_TEST_TARGETNAME} ${_PAR_ARGUMENTS})
  if( _TEST_BUILDS )

    # manual hacking includes
    if( DEFINED ${_TEST_NAME}_includedirs )
      include_directories(${${_TEST_NAME}_includedirs})
    endif()

    # separate exec
    add_executable( ${_TEST_TARGETNAME} ${_PAR_CPP} )

    # targets it depending on
    if( DEFINED _PAR_DEPENDS)
      add_dependencies( ${_TEST_TARGETNAME} ${_PAR_DEPENDS} )
    endif()


    # add dependency libraries if defined, proper way to do
    if( DEFINED _PAR_LIBS )
      target_link_libraries( ${_TEST_TARGETNAME} ${_PAR_LIBS} )
    endif()

    # add dependency libraries if defined, not so nice way
    if( DEFINED ${_TEST_NAME}_libraries )
      target_link_libraries(${${_TEST_NAME}_libraries})
    endif()

    # if mpi was found add it to the libraries
    if(_RUN_MPI)
      target_link_libraries( ${_TEST_TARGETNAME} ${MPI_LIBRARIES} )
    endif()

   IF (CF_ENABLE_UNITTESTS)
    # add boost unit test lib
    target_link_libraries( ${_TEST_TARGETNAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
   ENDIF()

    # prepare test
    if(_RUN_MPI)
      foreach(nprocs ${_PAR_MPI})
        if ( "${nprocs}" STREQUAL "default" )
          add_test(${_TEST_TARGETNAME}_dprocs ${CF_MPIRUN_PROGRAM} -np ${CF_TESTING_NB_PROCS} ${_TEST_COMMAND})
        else()
          add_test(${_TEST_TARGETNAME}_${nprocs}procs ${CF_MPIRUN_PROGRAM} -np ${nprocs} ${_TEST_COMMAND})
        endif()
      endforeach()
    else()
      add_test(  ${_TEST_TARGETNAME}_serial ${_TEST_COMMAND} )
    endif()

  endif( _TEST_BUILDS )

  # if installing
  if(CF_INSTALL_TESTS)
    install( TARGETS ${_TEST_TARGETNAME}
             RUNTIME DESTINATION ${CF_INSTALL_BIN_DIR}
             LIBRARY DESTINATION ${CF_INSTALL_LIB_DIR}
             ARCHIVE DESTINATION ${CF_INSTALL_LIB_DIR} )
  endif(CF_INSTALL_TESTS)

  # Set global variables
  set( ${_TEST_TARGETNAME}_builds  ${_TEST_BUILDS} )
  set( ${_TEST_TARGETNAME}_dir     ${_TEST_DIR} )
  set( ${_TEST_TARGETNAME}_libs    ${_PAR_LIBS} )
  if( DEFINED ${_TEST_TARGETNAME}_libs )
    list(REMOVE_DUPLICATES ${_TEST_TARGETNAME}_libs)
  endif()

  get_target_property( ${_TEST_TARGETNAME}_P_SOURCES   ${_TEST_TARGETNAME} SOURCES )
  get_target_property( ${_TEST_TARGETNAME}_LINK_FLAGS  ${_TEST_TARGETNAME} LINK_FLAGS )
  get_target_property( ${_TEST_TARGETNAME}_TYPE        ${_TEST_TARGETNAME} TYPE )

  # Log summary
#  MESSAGE(INFO "${_TEST_TARGETNAME}")
#  MESSAGE(INFO "    build requested     : [${CF_BUILD_${_TEST_TARGETNAME}}]")
#  MESSAGE(INFO "    builds              : [${${_TEST_TARGETNAME}_builds}]")
#  MESSAGE(INFO "    sources             : [${_PAR_CPP}]")
#  MESSAGE(INFO "    test dir            : [${_TEST_DIR}]")
#  MESSAGE(INFO "    mpirun              : [${CF_MPIRUN_PROGRAM}]")
#  MESSAGE(INFO "    mpi nb_proc         : [${_PAR_MPI}]")
#  MESSAGE(INFO "    libs                : [${${_TEST_TARGETNAME}_libs}]")
#  MESSAGE(INFO "    target type         : [${${_TEST_TARGETNAME}_TYPE}]")
#  MESSAGE(INFO "    target sources      : [${${_TEST_TARGETNAME}_P_SOURCES}]")
#  MESSAGE(INFO "    target link flags   : [${${_TEST_TARGETNAME}_LINK_FLAGS}]")

  # Add test to either CF_ENABLED_UTESTS  / CF_ENABLED_ATESTS
  #                 or CF_DISABLED_UTESTS / CF_DISABLED_ATESTS
  if(_TEST_BUILDS)
    if(_PAR_UTEST)
      set( CF_ENABLED_UTESTS ${CF_ENABLED_UTESTS} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    elseif(_PAR_PTEST)
      set( CF_ENABLED_PTESTS ${CF_ENABLED_PTESTS} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    endif()
  else()
    if(_PAR_UTEST)
      set( CF_DISABLED_UTESTS ${CF_DISABLED_UTESTS} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    elseif(_PAR_PTEST)
      set( CF_DISABLED_PTESTS ${CF_DISABLED_PTESTS} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    endif()
  endif()

#MESSAGE(FATAL_ERROR "BASTA")
endfunction( )

##############################################################################


