# try to compile using boost
IF ( NOT CF_HAVE_CUDA )
SET ( CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIR} )
CHECK_CXX_SOURCE_COMPILES (
	"#include <boost/math/special_functions/erf.hpp>
	int main(int argc, char* argv[])
	{
	 boost::math::erfc(0.);
	}"
	CF_HAVE_BOOST_ERFC )
ENDIF()

IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "MathTools" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND MathTools_files
FindMinimum.hh
MatrixInverter.cxx
MatrixEigenSolver.hh
IntersectSolver.hh
JacobiEigenSolver.cxx
LinearFunctor.hh
MathConsts.hh
MathFunctions.hh
ZeroDeterminantException.hh
ZeroDeterminantException.cxx
RealMatrix.hh
MatrixIntersect.hh
MatrixInverter.hh
InverterT.hh
OutOfBoundsException.hh
LUInverter.hh
RealVector.hh
ConstantFunctor.hh
InverterDiag.hh
LUInverter.cxx
MathChecks.hh
MathTools.hh
MatrixIntersect.cxx
MatrixEigenSolver.cxx
JacobiEigenSolver.hh
IntersectSolver.cxx
LUInverterT.hh
LUInverterT.ci
MatrixInverterT.hh
SVDInverter.hh
SVDInverter.cxx
RCM.h
RCM.cxx
CFMat.hh
CFVecSlice.hh
CFMatSlice.hh
ExprT.hh
MatExprT.hh
ArrayT.hh
MacrosET.hh
CFVec.hh
LeastSquaresSolver.cxx
LeastSquaresSolver.hh
# Function Parser (v4.5.2) from http://warp.povusers.org/FunctionParser/
FParser/fparser.cc
#FParser/fparser_gmpint.hh
FParser/fparser.hh
#FParser/fparser_mpfr.hh
FParser/fpconfig.hh
FParser/fpoptimizer.cc
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND MathTools_cflibs Common )
CF_ADD_KERNEL_LIBRARY ( MathTools )
ELSE()
FOREACH (AFILE ${MathTools_files} )
LIST(APPEND coolfluid-solver_files ../../src/MathTools/${AFILE} )
ENDFOREACH(AFILE)

CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()



