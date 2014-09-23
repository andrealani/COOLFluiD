#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"

#include "FiniteVolumeNavierStokes/NSFluxAndJacob.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolume/DerivativeComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSFluxAndJacob,
                       CellCenterFVMData,
		       ComputeDiffusiveFlux,
                       FiniteVolumeNavierStokesModule>
nsFluxAndJacobProvider("NSFluxAndJacobian");

//////////////////////////////////////////////////////////////////////////////

NSFluxAndJacob::NSFluxAndJacob(const std::string& name) :
  NSFlux<NavierStokesVarSet>(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NSFluxAndJacob::~NSFluxAndJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSFluxAndJacob::computeFlux(RealVector& result)
{ 
  if (_wallDistanceExists) {
      setWallDistance();
  }
  
  if (!this->getMethodData().isPerturb()) {
    GeometricEntity& geo = *getMethodData().getCurrentFace();
    SafePtr<DerivativeComputer> derivComputer = getMethodData().getDerivativeComputer();
    
    // this check is not needed
    // set the state values (pointers) corresponding to the
    // vertices of the control volume
    derivComputer->computeControlVolume(_states, &geo);
    _nbCVStates = derivComputer->getNbVerticesInControlVolume(&geo);
    
    _radius = 0.0;
    if (getMethodData().isAxisymmetric() || _isRadiusNeeded) {
      const Node& node0 = *geo.getNode(0);
      const Node& node1 = *geo.getNode(1);
      _radius = 0.5*(node0[YY] + node1[YY]);
    }
    
    // compute speed components and temperature in the given states
    // if you are updating in conservative variables your nodal values
    // MUST be already in primitive variables (there is inconsistency here !!!)
    _diffVar->setGradientVars(_states, _values, _nbCVStates);
    
    // compute control volume around the face and gradients
    derivComputer->computeGradients(&geo, _values, _gradients);
    
    // compute the average values
    derivComputer->computeAverageValues(&geo, _states, _avState);
    
    _diffVar->setComposition(_avState, false, this->getMethodData().iPerturbVar());
    
    const CFreal faceArea = this->socket_faceAreas.getDataHandle()[geo.getID()];
    
    // set the flux
    result = _diffVar->getFlux(_avState, _gradients, getMethodData().getUnitNormal(), _radius);
    result *= faceArea;
    
    // flux and its jacobian and flux are computed together
    const CFreal mu    = _diffVar->getCurrDynViscosity();
    const CFreal avRho = _diffVar->getDensity(_avState);
    const CFreal cvVolume = derivComputer->getControlVolume();
    const CFreal diffUpdateCoeff = mu*faceArea*faceArea/(avRho*cvVolume);

    DataHandle<CFreal> updateCoeff = this->socket_updateCoeff.getDataHandle();
    // left contribution to update coefficient
    const CFuint leftID = geo.getState(0)->getLocalID();
    updateCoeff[leftID] += diffUpdateCoeff;
    if (!geo.getState(1)->isGhost()) {
      // right contribution to update coefficient
      const CFuint rightID = geo.getState(1)->getLocalID();
      updateCoeff[rightID] += diffUpdateCoeff;
    }
    
    // jacobian computation
    SafePtr<vector<RealVector> > gradientsJacob = derivComputer->getGradientsJacob();
    RealVector& leftGradJacob  = (*gradientsJacob)[0];
    RealVector& rightGradJacob = (*gradientsJacob)[1];
    
    const RealVector& unitNormal = getMethodData().getUnitNormal();
    _diffVar->computeFluxJacobian(_avState, leftGradJacob, unitNormal, _radius, _lFluxJacobian);
    _lFluxJacobian *= faceArea; 
    
    if (!geo.getState(1)->isGhost()) {
      _diffVar->computeFluxJacobian(_avState, rightGradJacob, unitNormal, _radius, _rFluxJacobian);
      _rFluxJacobian *= faceArea;  
    }
  }
  else {
    // it can happen that a partial numerical perturbation is needed
    // in this case the base class has to be called
    NSFlux<NavierStokesVarSet>::computeFlux(result);
  }  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
