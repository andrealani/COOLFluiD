#include "FiniteElement/FiniteElement.hh"
#include "StdElementComputerStrategy.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/VolumeIntegrator.hh"

#include "ComputeConvectiveTerm.hh"
#include "ComputeDiffusiveTerm.hh"
#include "ComputeLinearSourceTerm.hh"
#include "ComputeIndepSourceTerm.hh"
#include "ComputeInertiaTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<StdElementComputerStrategy, FiniteElementMethodData, ComputeResidualStrategy, FiniteElementModule> stdElementComputerStrategyProvider("StdElementComputer");

//////////////////////////////////////////////////////////////////////////////

StdElementComputerStrategy::StdElementComputerStrategy(const std::string& name) :
  ComputeResidualStrategy(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StdElementComputerStrategy::~StdElementComputerStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdElementComputerStrategy::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ComputeResidualStrategy::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  _integResultMat.resize(nbEqs,nbEqs);
  _integResultVec.resize(nbEqs);

  // set up the Convective Term Computer
  SafePtr<ComputeConvectiveTerm> convTerm = getMethodData().getConvectiveTermComputer();
  if(!convTerm->isNull()) {
    CFLog(INFO,"Activating Convective Term: " << convTerm->getName() << "\n");
    _matrixComputeTerms.push_back(convTerm.d_castTo<ComputeTerm<FiniteElementMethodData> >());
  }

  // set up the Diffusive Term Computer
  SafePtr<ComputeDiffusiveTerm> diffTerm = getMethodData().getDiffusiveTermComputer();
  if(!diffTerm->isNull()) {
    CFLog(INFO,"Activating Diffusive Term: " << diffTerm->getName() << "\n");
    _matrixComputeTerms.push_back(diffTerm.d_castTo<ComputeTerm<FiniteElementMethodData> >());
  }
  // set up the Linear Source Term Computer
  SafePtr<ComputeLinearSourceTerm> linSTerm = getMethodData().getLinearSourceTermComputer();
  if(!linSTerm->isNull() && linSTerm->hasLinearCoef()) {
    CFLog(INFO,"Activating Linear Source Term: " << linSTerm->getName() << "\n");
    _matrixComputeTerms.push_back(linSTerm.d_castTo<ComputeTerm<FiniteElementMethodData> >());
  }

  // set up the Independent Source Term Computer
  SafePtr<ComputeIndepSourceTerm> indSTerm = getMethodData().getIndepSourceTermComputer();
  if(!indSTerm->isNull()) {
    CFLog(INFO,"Activating Independent Source Term: " << indSTerm->getName() << "\n");
    _vectorComputeTerms.push_back(indSTerm.d_castTo<ComputeTerm<FiniteElementMethodData> >());
  }

}

//////////////////////////////////////////////////////////////////////////////

void StdElementComputerStrategy::computeElemMatrix()
{
  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();

//   BlockAccumulator& acc = *local_elem_data.blockacc;
  RealMatrix& elemMat   = *local_elem_data.stiff_mat;
//   RealVector& elemVec   = *local_elem_data.load_vec;
  const CFuint nbEqs    =  local_elem_data.nbEqs;
  const CFuint nbStates =  local_elem_data.nbStates;
  GeometricEntity& cell = *local_elem_data.cell;

  vector<State*>& states = *cell.getStates();

  elemMat = 0.0;
/* CFout << "=======================================\n";
 CFout << "Cell: "<< cell.getID()<< "\n";
 CFout << "=======================================\n";*/
  vector<SafePtr<ComputeTerm<FiniteElementMethodData> > >::iterator term = _matrixComputeTerms.begin();
  for(; term != _matrixComputeTerms.end(); ++term)
  {
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {
      State& currState = *states[iState];
      local_elem_data.iState = iState;
      if (currState.isParUpdatable()) {
        for (CFuint jState = 0; jState < nbStates; ++jState) {
//           (*term)->setIndexes(iState,jState);
          local_elem_data.jState = jState;
          (*term)->computeTerm(&cell, _integResultMat);

          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
              elemMat(iState*nbEqs + iEq, jState*nbEqs + jEq) += _integResultMat(iEq,jEq);
            }
          }
        } // for jState
      } // for iState
    } // if isParUpdatable

  } // for term

}

//////////////////////////////////////////////////////////////////////////////

void StdElementComputerStrategy::computeElemVector()
{
  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();

//   BlockAccumulator& acc = *local_elem_data.blockacc;
//   RealMatrix& elemMat   = *local_elem_data.stiff_mat;
  RealVector& elemVec   = *local_elem_data.load_vec;
  const CFuint nbEqs    =  local_elem_data.nbEqs;
  const CFuint nbStates =  local_elem_data.nbStates;
  GeometricEntity& cell = *local_elem_data.cell;

  elemVec = 0.0;

  vector<SafePtr<ComputeTerm<FiniteElementMethodData> > >::iterator term = _vectorComputeTerms.begin();
  for(; term != _vectorComputeTerms.end(); ++term)
  {
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {

//       (*term)->setIndexes(iState,0);
      local_elem_data.iState = iState;
      local_elem_data.jState = 0;
      (*term)->computeTerm(&cell, _integResultVec);

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        elemVec[iState*nbEqs + iEq] += _integResultVec[iEq];
      }

    } // for iState
  } // for term
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

