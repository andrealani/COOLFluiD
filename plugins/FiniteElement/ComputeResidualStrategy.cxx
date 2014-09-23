#include "ComputeResidualStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

ComputeResidualStrategy::ComputeResidualStrategy(const std::string& name) : FiniteElementMethodStrategy(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeResidualStrategy::~ComputeResidualStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeResidualStrategy::computeElementResidual(std::vector<RealVector>& residual)
{

  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();

//   BlockAccumulator& acc = *local_elem_data.blockacc;
  RealMatrix& elemMat   = *local_elem_data.stiff_mat;
  RealVector& elemVec   = *local_elem_data.load_vec;
  const CFuint nbEqs    =  local_elem_data.nbEqs;
  const CFuint nbStates =  local_elem_data.nbStates;
  GeometricEntity& cell = *local_elem_data.cell;

  vector<State*>& states = *cell.getStates();

  cf_assert(elemMat.nbRows() == nbStates * nbEqs);
  cf_assert(elemMat.nbCols() == nbStates * nbEqs);
  cf_assert(elemVec.size()   == nbStates * nbEqs);
  cf_assert(residual.size()  <= nbStates);

  fill(residual.begin(),residual.end(),0.0);

  // compute element matrix and vector
  /// @todo consider joining all this functions together
  ///       into this one. Must then change ExplicitComputeSpaceResidual
  computeElemMatrix();
  computeElemVector();

  // compute the residual
  /// @todo try to make this product better
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      for (CFuint jState = 0; jState < nbStates; ++jState) {
        State& jStateRef = *(states[jState]);
        for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
          residual[iState][iEq] += elemMat(iState*nbEqs + iEq, jState*nbEqs + jEq) * jStateRef[jEq];
        }
      }
      residual[iState][iEq] -= elemVec[iState*nbEqs + iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeResidualStrategy::setup()
{
  CFAUTOTRACE;

  FiniteElementMethodStrategy::setup();

  // add specific configuration here
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
