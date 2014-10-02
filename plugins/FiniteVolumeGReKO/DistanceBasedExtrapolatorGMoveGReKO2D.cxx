#include "FiniteVolumeGReKO/DistanceBasedExtrapolatorGMoveGReKO2D.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "FiniteVolumeGReKO/FiniteVolumeGReKO.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveGReKO2D,
                       CellCenterFVMData,
                       NodalStatesExtrapolator<CellCenterFVMData>,
                       FiniteVolumeGReKOModule>
DistanceBasedExtrapolatorGMoveGReKO2DProvider("DistanceBasedGMoveGReKO2D");

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveGReKO2D::DistanceBasedExtrapolatorGMoveGReKO2D
(const std::string& name) :
  DistanceBasedExtrapolatorGMove(name),
  socket_wallDistance("wallDistance"),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveGReKO2D::~DistanceBasedExtrapolatorGMoveGReKO2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveGReKO2D::setup()
{
  DistanceBasedExtrapolatorGMove::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<NavierStokesTurb2DVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveGReKO2D::extrapolateInAllNodes()
{
  const CFuint OmegaID = 5;

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();

  const CFuint nbNodes = nodes.size();
  cf_assert (PhysicalModelStack::getActive()->getDim() == DIM_2D);

  const CFuint nbEqs =  PhysicalModelStack::getActive()->getNbEq();
  cf_assert (nbEqs == PhysicalModelStack::getActive()->getDim() + 6);


  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {

    const CFuint nbNeighborStates = _neighborStates[iNode].size();

    // reset to 0 the nodal states
    nodalStates[iNode] = 0.;
    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      const CFreal weight = static_cast<CFreal>(_weights[iNode][iState]);
      const State *const neighState = _neighborStates[iNode][iState];
      nodalStates[iNode] += (*neighState)*weight;
    }

    if (_isNodeToPrescribe[iNode]) {
      for (CFuint i = 0; i < _wValuesIdx.size(); ++i) {
        nodalStates[iNode][_wValuesIdx[i]] = _wValues[i];
      }

      //compute the reconstructed value for Omega
      // only internal States are considered
      nodalStates[iNode][OmegaID] = 0.;
      CFreal sumInnerWeight = 0.;
      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
        const State *const neighState = _neighborStates[iNode][iState];
        if (!neighState->isGhost()) {
          sumInnerWeight += static_cast<CFreal>(_weights[iNode][iState]);
        }
      }
      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
        const State *const neighState = _neighborStates[iNode][iState];

        if (!neighState->isGhost()) {
          //Compute distance to innerstate
          CFreal y0 = wallDistance[neighState->getLocalID()];

          //avoid too small distances
          y0 = max(y0, 10.e-10);

          const CFreal pdim =  (*neighState)[0] * _varSetTurb->getModel()->getPressRef();
          const CFreal Tdim =  (*neighState)[3] * _varSetTurb->getModel()->getTempRef();

          const CFreal mu = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
	    _diffVarTurb->getModel().getReferencePhysicalData()[NSTurbTerm::MU];

          const CFreal density = (*neighState)[0] / (_varSetTurb->getModel()->getR() * (*neighState)[3]);
          CFreal nu = mu / density;

          //this is not the best, but it avoids having to code another BC! because I
          //would have to dynamic cast to the GReKO varset to get the beta1
          const CFreal beta1 = 0.075;

          ///@todo here should this be adimensionalized (by the distance)???
          //Menter's definition
          const CFreal omegaWall = (10. * 6. * nu) / (beta1 * y0 * y0);

          const CFreal innerWeight = static_cast<CFreal>(_weights[iNode][iState])/sumInnerWeight;
          //omega are extrapolated only from inside the domain
          nodalStates[iNode][OmegaID] += omegaWall*innerWeight;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveGReKO2D::extrapolateInNodes
(const vector<Node*>& nodes)
{
  const CFuint OmegaID = 5;

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();

  const CFuint nbNodes = nodes.size();
  const CFuint nbEqs =  PhysicalModelStack::getActive()->getNbEq();
  cf_assert (nbEqs == PhysicalModelStack::getActive()->getDim() + 6);

  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const CFuint nodeID = nodes[iNode]->getLocalID();

    const CFuint nbNeighborStates = _neighborStates[nodeID].size();

    // reset to 0 the nodal states
    nodalStates[nodeID] = 0.;
    for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
      const CFreal weight = static_cast<CFreal>(_weights[nodeID][iState]);
      const State *const neighState = _neighborStates[nodeID][iState];
      nodalStates[nodeID] += (*neighState)*weight;
    }

    if (_isNodeToPrescribe[nodeID]) {
      for (CFuint i = 0; i < _wValuesIdx.size(); ++i) {
        nodalStates[nodeID][_wValuesIdx[i]] = _wValues[i];
      }

      //compute the reconstructed value for Omega
      // only internal States are considered
      nodalStates[nodeID][OmegaID] = 0.;
      CFreal sumInnerWeight = 0.;
      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
        const State *const neighState = _neighborStates[nodeID][iState];
        if (!neighState->isGhost()) {
          sumInnerWeight += static_cast<CFreal>(_weights[nodeID][iState]);
        }
      }

      for (CFuint iState = 0; iState < nbNeighborStates; ++iState) {
        const State *const neighState = _neighborStates[nodeID][iState];
        if (!neighState->isGhost()) {
          //Compute distance to innerstate
          CFreal y0 = wallDistance[neighState->getLocalID()];

          //avoid too small distances
          y0 = max(y0, 10.e-10);

          const CFreal pdim =  (*neighState)[0] * _varSetTurb->getModel()->getPressRef();
          const CFreal Tdim =  (*neighState)[3] * _varSetTurb->getModel()->getTempRef();

          const CFreal mu = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
	    _diffVarTurb->getModel().getReferencePhysicalData()[NSTurbTerm::MU];

          const CFreal density = (*neighState)[0] / (_varSetTurb->getModel()->getR() * (*neighState)[3]);
          CFreal nu = mu / density;

          //this is not the best, but it avoids having to code another BC! because I
          //would have to dynamic cast to the GReKO varset to get the beta1
          const CFreal beta1 = 0.075;

          ///@todo here should this be adimensionalized (by the distance)???
          //Menter's definition
          const CFreal omegaWall = (10. * 6. * nu) / (beta1 * y0 * y0);

          const CFreal innerWeight = static_cast<CFreal>(_weights[iNode][iState])/sumInnerWeight;
          //omega are extrapolated only from inside the domain
          nodalStates[nodeID][OmegaID] += omegaWall*innerWeight;
        }
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
