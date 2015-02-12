# Function to add a test.
#
# Mandatory keywords:
# - UCASE/PCASE
#      choose according to the test profile: unit-test, performance-test
#      give as value the relative_to_source dir/name of the CFcase
#
# Optional keywords:
# - ARGUMENTS
#      a list of arguments to pass to the the test execution (default empty)
# - MPI (default: if HAVE_MPI then runs the case on default number of processors)
#      number of processors to use for call with mpirun (it can be a list of numbers)
#      the keyword "default" will take the global number CF_TESTING_NB_PROCS
#
# The two master switches turns each type on and off
#   - CF_ENABLE_UNIT_CASES
#   - CF_ENABLE_PERFORMANCE_CASES
#
# After calling this function, the test is added to one of the following lists:
#   - CF_ENABLED_UCASES
#   - CF_DISABLED_UCASES
#   - CF_ENABLED_PCASES
#   - CF_DISABLED_PCASES

function( cf_add_case )

#MESSAGE(INFO,"***CASE*************************************************************")
#MESSAGE(INFO,"***CASE*** CF_ENABLE_UNIT_CASES= ${CF_ENABLE_UNIT_CASES}")
#MESSAGE(INFO,"***CASE*** CF_ENABLE_PERFORMANCE_CASES= ${CF_ENABLE_PERFORMANCE_CASES}")

  set( single_value_args UCASE PCASE CASEDIR )
#  set( multi_value_args  MPI)
  set( multi_value_args  MPI CASEFILES)
  
  set( _TEST_DIR ${CMAKE_CURRENT_BINARY_DIR} )

#MESSAGE(INFO,"***CASE*** CF_ENABLE_PERFORMANCE_CASES= ${CF_ENABLE_PERFORMANCE_CASES}")

  # parse and complain if stg wrong with the arguments
  cmake_parse_arguments(_PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN})

  SET( _BUILD_CASEDIR  ${CMAKE_CURRENT_BINARY_DIR}/${_PAR_CASEDIR} )
 
  # create testcase folder  
  FILE(MAKE_DIRECTORY [${_BUILD_CASEDIR}])

  # copy all needed files into the testcase folder   
  FOREACH( ACFG ${_PAR_CASEFILES} )
   FILE ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${_PAR_CASEDIR}/${ACFG} DESTINATION ${_BUILD_CASEDIR}) 
   LOGVERBOSE("***Copying ${CMAKE_CURRENT_SOURCE_DIR}/${_PAR_CASEDIR}/${ACFG} to ${_BUILD_CASEDIR}")
  ENDFOREACH()

  foreach(unparsed_arg ${_PAR_UNPARSED_ARGUMENTS})
    if( NOT ( (${unparsed_arg} STREQUAL "AND") OR (${unparsed_arg} STREQUAL "OR") OR (${unparsed_arg} STREQUAL "NOT") ) )
      list( APPEND _UNPARSED_ARGUMENTS ${unparsed_arg} )
    endif()
  endforeach()
  if(_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unknown keywords given to cf_add_test(): \"${_UNPARSED_ARGUMENTS}\"")
  endif()
  if( (NOT _PAR_UCASE) AND (NOT _PAR_PCASE))
    message(FATAL_ERROR "The call to cf_add_test() doesn't set the required \"UCASE/ATEST/PCASE test-name\" argument.")
  endif()

  # select unit/acceptance/performance tests, removing extension if present
  if(_PAR_UCASE)
    set(_TEST_NAME ${_PAR_CASEDIR}/${_PAR_UCASE})
    set(_TEST_TARGETNAME "case-unit-")
    set(_TEST_ENABLED ${CF_ENABLE_UNIT_CASES})
  elseif(_PAR_PCASE)
    set(_TEST_NAME ${_PAR_CASEDIR}/${_PAR_PCASE})
    set(_TEST_TARGETNAME "case-perf-")
    set(_TEST_ENABLED ${CF_ENABLE_PERFORMANCE_CASES})
  endif()
  string( REPLACE ".CFcase"  "" _TEST_NAME ${_TEST_NAME} )
  string( REPLACE ".XCFcase" "" _TEST_NAME ${_TEST_NAME} )

  # enrich test name with path
  set(_TEST_WDIR ${CMAKE_CURRENT_BINARY_DIR}/${_TEST_NAME})
  set(_TEST_NAME ${CMAKE_CURRENT_SOURCE_DIR}/${_TEST_NAME})
  set(_TEST_CFCASE ${_TEST_NAME}.CFcase)
  set(_TEST_CFCASE_SHORT ${_TEST_CFCASE})
  string(REPLACE "${CMAKE_SOURCE_DIR}/" "" _TEST_CFCASE_SHORT ${_TEST_CFCASE_SHORT})
  string(REPLACE "${CMAKE_SOURCE_DIR}/" "" _TEST_NAME ${_TEST_NAME})
  string(REPLACE "/" "-" _TEST_NAME ${_TEST_NAME})
  set(_TEST_TARGETNAME "${_TEST_TARGETNAME}${_TEST_NAME}")


#MESSAGE(INFO,"***CASE*** _TEST_NAME= ${_TEST_NAME}")
#MESSAGE(INFO,"***CASE*** _TEST_CFCASE= ${_TEST_CFCASE}")
#MESSAGE(INFO,"***CASE*** _TEST_TARGETNAME= ${_TEST_TARGETNAME}")
#MESSAGE(INFO,"***CASE*** _CMAKE_CURRENT_SOURCE_DIR= ${CMAKE_CURRENT_SOURCE_DIR}")
#MESSAGE(INFO,"***CASE*** _CMAKE_SOURCE_DIR= ${CMAKE_SOURCE_DIR}")
#MESSAGE(INFO,"***CASE*** _TEST_ENABLED= ${_TEST_ENABLED}")

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

  # build the testing command
  set(_TEST_COMMAND "${coolfluid_solver_exe}" "--scase" "${CMAKE_BINARY_DIR}/${_TEST_CFCASE_SHORT}")
  if( _TEST_BUILDS )

    # read the file and chew it into a nice formatting
    file(READ "${_TEST_CFCASE}" contents)
    string(      REPLACE  ";" "__REPLACED_SEMICOLON__" contents "${contents}") # semicolons are the dividers between list items
    string(REGEX REPLACE  "[\\][ \t]*[\n]" " " contents "${contents}") # multiline text broken by backslash
    string(REGEX REPLACE  "[ \t]*[\n]+[ \t]*" "\n" contents "${contents}") # clearing whitespaces at the end of line
    string(REGEX REPLACE  "^[ \t]+" "" contents "${contents}") # previous skips possible whitespaces in first line
    string(REGEX REPLACE  "[ \t]+([=]*)[ \t]*" "\\1 " contents "${contents}") # everywhere always exactly one whitespace, plus clearing in front of = sign
    # preparing the list containing the blocks, 
    # the keyword is #!<space>KEYWORD, rest of the line must be ignored (until first newline)
    # replace back any __REPLACED_SEMICOLON__ symbols at the end
    string(      REPLACE  "#!" ";" blockcontents "${contents}") # replacing terminating #! to ; (this turns it into a list)
    string(REGEX MATCHALL "[!][#][^;]*[;]" blockcontents "${blockcontents}") # keeping !# ... ; enclosures
    # further formatting
    foreach(iblock ${blockcontents})
      string(    REPLACE  "${iblock}" " " contents "${contents}") # removing blocks
    endforeach()
    string(REGEX REPLACE "[\n]" ";" contents "${contents}") # now each line is an item in the list
    # exporting lines having ###, first wor is key, 
    # equation sign is attached (if present)
    string(REGEX MATCHALL "[#][#][#][^;]*[;]" pragmacontents "${contents}") # now each line is an item in the list
    string(REGEX REPLACE  "[#][#][#][ \t]*" "" pragmacontents "${pragmacontents}") # now each line is an item in the list
    # pure cfcase options
    string(REGEX REPLACE "[ \t]*[#][^;]*[;]+" ";" settingcontents "${contents}") # clear everything in a line after any hash, including the hash
    string(REGEX REPLACE "[=][ \t]*" "=" settingcontents "${settingcontents}") # just for convenience, then one can split by space

    # process pragmas
    foreach(ipragma ${pragmacontents})

      # Residual check
      if ( ipragma MATCHES "Residual=.*" ) 
       string(REGEX REPLACE "Residual=[ \t]*" "" CFtargetResidual ${ipragma}) 
       string( REGEX REPLACE "[-+0-9eE. \t]" "" CFcheckresidual "${CFtargetResidual}")
        if( NOT ( "${CFcheckresidual}" STREQUAL "" ) )
           message(FATAL_ERROR "Syntax error for testcase residual: ${_TEST_CFCASE_SHORT}")
        endif()
        list(APPEND _TEST_COMMAND "--residual" "${CFtargetResidual}")
      endif()

    endforeach(ipragma ${pragmacontents})

#MESSAGE(INFO "***CASE*** RESIDUAL=${CFtargetResidual}")

    # process cfoptions
    foreach(icfopt ${settingcontents})

      # extract working directory
      if( icfopt MATCHES "Simulator.Paths.WorkingDir=.*" )
	string( REGEX REPLACE "Simulator.Paths.WorkingDir=(.*)" "\\1" CFtestWorkdir "${icfopt}")
	string(SUBSTRING "${CFtestWorkdir}" 0 1 test_local)
	if ( "${test_local}" STREQUAL ".")
          message(FATAL_ERROR "Testcase's working directory is local (${_TEST_CFCASE_SHORT}), please adjust Simulator.Paths.WorkingDir in the CFcase.")
	endif()
	file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/${CFtestWorkdir}")
      endif()

      # extract libraries and check if compiled
      if( icfopt MATCHES "Simulator.Modules.Libs=.*" )
	STRING ( REGEX REPLACE "Simulator.Modules.Libs=(.*)" "\\1" CFtestLibs "${icfopt}")
        separate_arguments(CFtestLibs)
	set(CFmissingLibs "")
	foreach( ilib ${CFtestLibs} )
	  string( REGEX REPLACE ".*lib([-a-zA-Z0-9]+)" "\\1" ilibname "${ilib}")
         #message (WARNING "Testing library message ${ilibname}"  ) 
	  if( NOT CF_COMPILES_${ilibname} )
            set(CFmissingLibs "${CFmissingLibs} ${ilibname}")
	  endif()
	endforeach()
	if( NOT ( "${CFmissingLibs}" STREQUAL "" ) )
	  message(WARNING "Test case is turned off (${_TEST_CFCASE_SHORT}) because the following libs are missing: ${CFmissingLibs}")
	  set(_TEST_BUILDS FALSE)
	endif()
      endif()

    endforeach(icfopt ${settingcontents})

    # process blocks like latex description or file list
    #foreach(iblock ${blockcontents})
    #endforeach(iblock ${blockcontents})

#MESSAGE(INFO "***CASE*** WORKDIR=|${CFtestWorkdir}|")

    # additional options
    list(APPEND _TEST_COMMAND "--tolerance" "5" )
    list(APPEND _TEST_COMMAND "--bdir" "${COOLFluiD_SOURCE_DIR}" )
#    list(APPEND _TEST_COMMAND "--bdir" "${COOLFluiD_BINARY_DIR}" )
    list(APPEND _TEST_COMMAND "--ldir" "${COOLFluiD_BINARY_DIR}/dso" )

  endif( _TEST_BUILDS )

#MESSAGE("TESTCOMMAND  ${_TEST_COMMAND}")

  # prepare
  if( _TEST_BUILDS )
  
    # time to copy the CFcase to binary tree
    CONFIGURE_FILE ( ${CMAKE_SOURCE_DIR}/${_TEST_CFCASE_SHORT} ${CMAKE_BINARY_DIR}/${_TEST_CFCASE_SHORT} @ONLY )

    # prepare test
    if(_RUN_MPI)
      foreach(nprocs ${_PAR_MPI})
        if ( "${nprocs}" STREQUAL "default" )
          add_test(NAME ${_TEST_TARGETNAME}_dprocs COMMAND ${CF_MPIRUN_PROGRAM} -np ${CF_TESTING_NB_PROCS} ${_TEST_COMMAND})
          #add_test(NAME ${_TEST_TARGETNAME}_dprocs WORKING_DIRECTORY ${_TEST_WDIR} COMMAND ${CF_MPIRUN_PROGRAM} -np ${CF_TESTING_NB_PROCS} ${_TEST_COMMAND})
        else()
          add_test(NAME ${_TEST_TARGETNAME}_${nprocs}procs COMMAND ${CF_MPIRUN_PROGRAM} -np ${nprocs} ${_TEST_COMMAND})
          #add_test(NAME ${_TEST_TARGETNAME}_${nprocs}procs WORKING_DIRECTORY ${_TEST_WDIR} COMMAND ${CF_MPIRUN_PROGRAM} -np ${nprocs} ${_TEST_COMMAND})
        endif()
      endforeach()
    else()
      add_test(NAME   ${_TEST_TARGETNAME}_serial COMMAND ${_TEST_COMMAND})
      #add_test(NAME   ${_TEST_TARGETNAME}_serial WORKING_DIRECTORY ${_TEST_WDIR} COMMAND ${_TEST_COMMAND})
   endif()

  endif( _TEST_BUILDS )

  # if installing
  if(CF_INSTALL_TESTS)
    message(WARNING "Installing Testcases are not implemented.")
  endif(CF_INSTALL_TESTS)

  # Add test to either CF_ENABLED_UCASES  / CF_ENABLED_PCASES
  #                 or CF_DISABLED_UCASES / CF_DISABLED_PCASES
  if(_TEST_BUILDS)
    if(_PAR_UCASE)
      set( CF_ENABLED_UCASES ${CF_ENABLED_UCASES} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    elseif(_PAR_PCASE)
      set( CF_ENABLED_PCASES ${CF_ENABLED_PCASES} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    endif()
  else()
    if(_PAR_UCASE)
      set( CF_DISABLED_UCASES ${CF_DISABLED_UCASES} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    elseif(_PAR_PCASE)
      set( CF_DISABLED_PCASES ${CF_DISABLED_PCASES} ${_TEST_TARGETNAME} CACHE INTERNAL "" )
    endif()
  endif()

#MESSAGE(FATAL_ERROR "BASTA")
endfunction( )

##############################################################################


