#include "Framework/PhysicalModel.hh"
#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/LUFactorization.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LUFactorization, LUSGSIteratorData, LUSGSMethodModule> luFactorizationProvider("LUFact");


//////////////////////////////////////////////////////////////////////////////

LUFactorization::LUFactorization(std::string name) :
  LUSGSIteratorCom(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  socket_isStatesSetParUpdatable("isStatesSetParUpdatable"),
  socket_statesSetIdx("statesSetIdx"),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

void LUFactorization::execute()
{
  CFAUTOTRACE;

  // Gets the diagonal matrices
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();
  cf_assert(diagBlockJacobMatr.size() == getMethodData().getNbrStatesSets());

  // Gets the number of states sets
  const CFuint nbrStatesSets = diagBlockJacobMatr.size();

  // get isStatesSetParUpdatable data handle
  DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();

  // Loops over the states sets
  for (CFuint iSet = 0; iSet < nbrStatesSets; ++iSet)
  {
    if (isStatesSetParUpdatable[iSet])
    {
      // Dereferences the current matrix
      RealMatrix& currDiagMatrix = diagBlockJacobMatr[iSet];

      // number of states in the current set
      const CFuint resSize = currDiagMatrix.nbRows();
      
      // actual LU factorization
      // loop over the diagonal elements
      const CFuint resSizeM1 = resSize - 1;
      for (CFuint iDiag = 0; iDiag < resSizeM1; ++iDiag)
      {
        // LU FACTORIZATION
        factorizeMatrix(iDiag,currDiagMatrix);
      }
    }
  }

  // set second element of socket_statesSetIdx to 0 --> updateCoefs are not recomputed
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  statesSetIdx[1] = 0;
}

//////////////////////////////////////////////////////////////////////////////

void LUFactorization::factorizeMatrix(const CFuint diag, RealMatrix& matrix)
{
  const CFuint size = matrix.nbRows();
  const CFreal invDiag = 1.0/matrix(diag,diag);
  for (CFuint iRow = diag+1; iRow < size; ++iRow)
  {
    const CFreal factor = matrix(iRow,diag)*invDiag;
    matrix(iRow,diag) = factor;// L matrix
    for (CFuint iCol = diag+1; iCol < size; ++iCol)
    {
      matrix(iRow,iCol) -= factor*matrix(diag,iCol);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LUFactorization::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  LUSGSIteratorCom::setup();

  // get the number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > LUFactorization::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_diagBlockJacobMatr     );
  result.push_back(&socket_isStatesSetParUpdatable);
  result.push_back(&socket_statesSetIdx           );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
