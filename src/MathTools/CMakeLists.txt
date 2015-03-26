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
FunctionParser.hh
FunctionParser.cxx
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
LeastSquaresSolver.hh
)

LIST ( APPEND MathTools_cflibs Common )

CF_ADD_KERNEL_LIBRARY ( MathTools )

CF_WARN_ORPHAN_FILES()



