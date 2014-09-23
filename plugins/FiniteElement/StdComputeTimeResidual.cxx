#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/BlockAccumulator.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/StdComputeTimeResidual.hh"
#include "FiniteElement/ComputeInertiaTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdComputeTimeResidual, FiniteElementMethodData, FiniteElementModule> stdComputeTimeResidualProvider("StdComputeTimeResCom");

//////////////////////////////////////////////////////////////////////////////

StdComputeTimeResidual::StdComputeTimeResidual(const std::string& name) :
  FiniteElementMethodCom(name),
  socket_updateCoeff("updateCoeff"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

StdComputeTimeResidual::~StdComputeTimeResidual()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdComputeTimeResidual::setup()
{
  // first call parent method
  FiniteElementMethodCom::setup();

}

//////////////////////////////////////////////////////////////////////////////

void StdComputeTimeResidual::executeOnTrs()
{
  CFAUTOTRACE;

/// @todo Relaxation: fix it or leave it
#if 0 // experience with relaxation

  SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLinearSystemSolver()->getMatrix();
  // get the data handle to the update coeff
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();

  updateCoeff = 0.0;

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  for(CFuint iTR = 0; iTR < trs->getNbTRs(); ++iTR) {

    TopologicalRegion * tr = trs->getTopologicalRegion(iTR);

    for(CFuint iGeoEnt = 0; iGeoEnt < tr->getNbGeomEnt(); ++iGeoEnt) {

      CFLogDebugMax("Cell " << iGeoEnt << "\n");

      GeometricEntity& cell = *tr->getGeomEnt(iGeoEnt);
      vector<State*>& cellStates = *cell.getStates();
      const CFuint nbStatesInCell = cellStates.size();

      const CFreal nodeVolume = cell.computeVolume() / nbStatesInCell;

      CFLogDebugMax("Cell Volume" << nodeVolume*nbStatesInCell << "\n");

      // write to the right hand side
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        if(cellStates[iState]->isParUpdatable()) {
          const CFuint stateID = cellStates[iState]->getLocalID();
            updateCoeff[stateID] += nodeVolume;
        }
      }
    }
  }

  LSSIdxMapping& idxMapping = LSSIdxMapping::getInstance();

  // add the diagonal entries in the jacobian
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    if(states[iState]->isParUpdatable()) {

      const CFuint stateID = states[iState]->getLocalID();
      CFuint globalID = idxMapping.getID(stateID)*nbEqs;

      const CFreal diagValue = updateCoeff[stateID]/cfl;

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq, ++globalID) {
        jacobMatrix->addValue(globalID, globalID, diagValue);
      }
    }
  }
#endif
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdComputeTimeResidual::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
