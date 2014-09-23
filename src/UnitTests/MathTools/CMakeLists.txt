  LIST ( APPEND TestSuite_MathTools_libs MathTools)
  LIST ( APPEND TestSuite_MathTools_files
    TestSuite_MathTools.cxx
    Test_MatrixInverter.hh
  )

#  CF_ADD_TEST ( TestSuite_MathTools )
  
  # TODO: Test_RCuthillMckee unit test must be merged togetherwitht the testsuite
  LIST ( APPEND Test_RCuthillMckee_libs   MathTools )
  LIST ( APPEND Test_RCuthillMckee_files  Test_RCuthillMckee.cpp )
  
#  CF_ADD_TEST ( Test_RCuthillMckee )

  CONFIGURE_FILE ( struct_1x10.txt ${CMAKE_CURRENT_BINARY_DIR}/struct_1x10.txt @ONLY )
  CONFIGURE_FILE ( struct_3x10.txt ${CMAKE_CURRENT_BINARY_DIR}/struct_3x10.txt @ONLY )
  CONFIGURE_FILE ( testmatrix.txt  ${CMAKE_CURRENT_BINARY_DIR}/testmatrix.txt  @ONLY )
  CONFIGURE_FILE ( testmesh.txt    ${CMAKE_CURRENT_BINARY_DIR}/testmesh.txt    @ONLY )


cf_add_test(
  UTEST mathtools1
  CPP   TestSuite_MathTools.cxx
        Test_MatrixInverter.hh
  LIBS  MathTools
)


cf_add_test(
  UTEST mathtools2
  CPP   TestSuite_MathTools.cxx
        Test_MatrixInverter.hh
  LIBS  MathTools
  MPI   2 default
)
cf_add_test(
  PTEST mathtools3
  CPP   TestSuite_MathTools.cxx
        Test_MatrixInverter.hh
  LIBS  MathTools
  MPI   1 4
)

#cf_add_testcase (
#  UCASE name
#  CFCASE *****.CFcase torch.inter
#  MPI   1


  CF_WARN_ORPHAN_FILES()

