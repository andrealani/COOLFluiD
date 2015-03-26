LIST ( APPEND TestSuite_MathTools_libs MathTools)

LIST ( APPEND TestSuite_MathTools_files
utest-leastSquaresSolver.cxx  
utest-matrixInverter.cxx	
utest-realVector.cxx
)

cf_add_test(
  UTEST matrixInverter
  CPP   utest-matrixInverter.cxx 
  LIBS  MathTools
)


cf_add_test(
  UTEST realVector 
  CPP   utest-realVector.cxx
  LIBS  MathTools
)

cf_add_test(
  UTEST leastSquaresSolver
  CPP   utest-leastSquaresSolver.cxx
  LIBS  MathTools
)

CF_WARN_ORPHAN_FILES()

