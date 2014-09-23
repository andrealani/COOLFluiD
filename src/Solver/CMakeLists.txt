### coolfluid-solver ##########################################################

CONFIGURE_FILE ( coolfluid-solver.xml ${CMAKE_CURRENT_BINARY_DIR}/coolfluid-solver.xml       @ONLY )
CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/coolfluid-solver-wrapper ${CMAKE_CURRENT_BINARY_DIR}/coolfluid-solver-wrapper @ONLY )

LIST ( APPEND coolfluid-solver_files coolfluid-solver.cxx )
LIST ( APPEND coolfluid-solver_libs ${CF_KERNEL_LIBS} ${CF_KERNEL_STATIC_LIBS} ${CF_Boost_LIBRARIES} )

CF_ADD_PLUGIN_APP ( coolfluid-solver )

# testing of testcases
GET_PROPERTY ( coolfluid_solver_path TARGET coolfluid-solver PROPERTY LOCATION )
SET ( coolfluid_solver_exe ${coolfluid_solver_path} CACHE "Full path to coolfluid-solver" INTERNAL )
MARK_AS_ADVANCED ( coolfluid_solver_exe )

###############################################################################


CF_WARN_ORPHAN_FILES()

