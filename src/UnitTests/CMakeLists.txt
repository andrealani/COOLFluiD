# little helping executable for comparing CFmeshes
LIST ( APPEND test-tools-cfmesh-compare_files test-tools-cfmesh-compare.cxx )
LIST ( APPEND test-tools-cfmesh-compare_libs ${CF_KERNEL_LIBS} ${CF_KERNEL_STATIC_LIBS} )

CF_ADD_PLUGIN_APP ( test-tools-cfmesh-compare )

GET_PROPERTY ( test-tools-cfmesh-compare_path TARGET test-tools-cfmesh-compare PROPERTY LOCATION )
SET ( test-tools-cfmesh-compare_exe ${test-tools-cfmesh-compare_path} CACHE "Full path to test-tools-cfmesh-compare" INTERNAL )
MARK_AS_ADVANCED ( test-tools-cfmesh-compare_exe )

IF (NOT CF_HAVE_CUDA)
add_subdirectory ( MathTools )
ENDIF()
