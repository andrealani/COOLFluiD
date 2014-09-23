#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/ComputeStatesSetUpdatePivot.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeStatesSetUpdatePivot, LUSGSIteratorData, LUSGSMethodModule>
    computeStatesSetUpdatePivotProvider("ComputeStatesSetUpdatePivot");


//////////////////////////////////////////////////////////////////////////////

ComputeStatesSetUpdatePivot::ComputeStatesSetUpdatePivot(std::string name) :
  ComputeStatesSetUpdate(name),
  socket_pivotLUFactorization("pivotLUFactorization"),
  m_resAux(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeStatesSetUpdatePivot::execute()
{
  // Gets current states set index
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  const CFuint currIdx = statesSetIdx[0];

  // get isStatesSetParUpdatable data handle
  DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();

  if (isStatesSetParUpdatable[currIdx])
  {
    // Gets the diagonal matrices
    DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

    // Gets the pivoting vectors
    DataHandle< vector< CFuint > > pivotLUFactorization = socket_pivotLUFactorization.getDataHandle();

    // Gets the rhs vectors
    DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

    // Dereferences the current matrix
    RealMatrix& currDiagMatrix = diagBlockJacobMatr[currIdx];

    //Dereferences the pivot elements
    vector< CFuint >& currPivot = pivotLUFactorization[currIdx];

    // Compute the size of the residuals
    const CFuint currResSize = currPivot.size();
    cf_assert(rhsCurrStatesSet.size() >= currResSize);
    cf_assert(rhsCurrStatesSet.size() >= currDiagMatrix.nbRows());
    cf_assert(currResSize == currDiagMatrix.nbRows());

    // Derefence the auxiliary rhs variable
    RealVector& resAux = *m_resAux;
    cf_assert(rhsCurrStatesSet.size() == resAux.size());

    // Rearrange the elements of the rhsCurrStatesSet vector. resAux is used to hold them.
    // At the end of the algorithm rhsCurrStatesSet vector will contain the current states set update.
    for (CFuint iRes = 0; iRes < currResSize; ++iRes)
    {
      const CFuint jRes = currPivot[iRes];
      resAux[iRes] = rhsCurrStatesSet[jRes];
    }

    // Solve the two triangular systems
    solveTriangularSystems(currDiagMatrix, resAux);

    // Copy resAux into rhsCurrStatesSet
    for (CFuint iRes = 0; iRes < currResSize; ++iRes)
    {
      rhsCurrStatesSet[iRes] = resAux[iRes];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeStatesSetUpdatePivot::setup()
{
  CFAUTOTRACE;

  // get the auxiliary rhs variable
  m_resAux = getMethodData().getResAux();

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ComputeStatesSetUpdatePivot::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ComputeStatesSetUpdate::needsSockets();

  result.push_back(&socket_pivotLUFactorization);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
