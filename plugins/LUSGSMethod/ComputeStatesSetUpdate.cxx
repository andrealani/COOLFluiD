#include "LUSGSMethod/LUSGSMethod.hh"
#include "LUSGSMethod/ComputeStatesSetUpdate.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeStatesSetUpdate, LUSGSIteratorData, LUSGSMethodModule>
    computeStatesSetUpdateProvider("ComputeStatesSetUpdate");


//////////////////////////////////////////////////////////////////////////////

ComputeStatesSetUpdate::ComputeStatesSetUpdate(std::string name) :
  LUSGSIteratorCom(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_statesSetIdx("statesSetIdx"),
  socket_isStatesSetParUpdatable("isStatesSetParUpdatable"),
  m_resAux()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeStatesSetUpdate::execute()
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

    // Gets the rhs vectors
    DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();

    // Dereferences the current matrix
    RealMatrix& currDiagMatrix = diagBlockJacobMatr[currIdx];


    // Compute the size of the residuals
    const CFuint currResSize = currDiagMatrix.nbRows();
    cf_assert(rhsCurrStatesSet.size() >= currResSize);
    cf_assert(rhsCurrStatesSet.size() >= currDiagMatrix.nbRows());

    // Derefence the auxiliary rhs variable
    RealVector& resAux = *m_resAux;
    cf_assert(rhsCurrStatesSet.size() == resAux.size());

    // Copy rhsCurrStatesSet into resAux
    for (CFuint iRes = 0; iRes < currResSize; ++iRes)
    {
      resAux[iRes] = rhsCurrStatesSet[iRes];
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

void ComputeStatesSetUpdate::solveTriangularSystems(const RealMatrix& lhsMatrix, RealVector& rhs)
{
  const CFuint size = lhsMatrix.nbRows();
  cf_assert(size <= rhs.size());

  // Solve the lower triangular system
  // Do forward substitution
  for (CFuint i = 1; i < size; ++i)
  {
    for (CFuint j = 0; j < i; ++j)
    {
      rhs[i] -= lhsMatrix(i,j)*rhs[j];
    }
  }

  // Solve the upper triangular system
  // Do backward substitution
  const CFuint sizeM1 = size-1;
  rhs[sizeM1] /= lhsMatrix(sizeM1,sizeM1);
  for (CFint i = size - 2; i >= 0; --i)
  {
    for (CFuint j = i + 1; j < size; ++j)
    {
      rhs[i] -= lhsMatrix(i,j)*rhs[j];
    }
    rhs[i] /= lhsMatrix(i,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeStatesSetUpdate::setup()
{
  CFAUTOTRACE;

  // get the auxiliary rhs variable
  m_resAux = getMethodData().getResAux();

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ComputeStatesSetUpdate::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_diagBlockJacobMatr);
  result.push_back(&socket_rhsCurrStatesSet);
  result.push_back(&socket_statesSetIdx);
  result.push_back(&socket_isStatesSetParUpdatable);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD
