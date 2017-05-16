### coolfluid-solver ##########################################################

CONFIGURE_FILE ( coolfluid-solver.xml ${CMAKE_CURRENT_BINARY_DIR}/coolfluid-solver.xml       @ONLY )
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/coolfluid-solver-wrapper ${CMAKE_CURRENT_BINARY_DIR}/coolfluid-solver-wrapper @ONLY )

LIST ( APPEND coolfluid-solver_files coolfluid-solver.cxx )

LOG("coolfluid-solver_files is: ${coolfluid-solver_files}")
LIST ( APPEND coolfluid-solver_libs ${CF_KERNEL_LIBS} ${CF_KERNEL_STATIC_LIBS} )
 
IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND coolfluid-solver_libs ${CF_Boost_LIBRARIES} )
ELSE()
LIST ( APPEND coolfluid-solver_libs dl )
ENDIF()

MESSAGE ("###### coolfluid-solver_libs are ${coolfluid-solver_libs}")
MESSAGE ("###### coolfluid-solver_cflibs are ${coolfluid-solver_cflibs}")
IF ( CF_HAVE_ALLSTATIC )
LOG ("coolfluid-solver_libs are ${coolfluid-solver_libs}")
LOG("CF_KERNEL_LIBS is [${CF_KERNEL_LIBS}]" )
LOG("CF_KERNEL_STATIC_LIBS is [${CF_KERNEL_STATIC_LIBS}]" )
ENDIF()

CF_ADD_PLUGIN_APP ( coolfluid-solver )

# testing of testcases
GET_PROPERTY ( coolfluid_solver_path TARGET coolfluid-solver PROPERTY LOCATION )
SET ( coolfluid_solver_exe ${coolfluid_solver_path} CACHE "Full path to coolfluid-solver" INTERNAL )
MARK_AS_ADVANCED ( coolfluid_solver_exe )

###############################################################################

CF_WARN_ORPHAN_FILES()

