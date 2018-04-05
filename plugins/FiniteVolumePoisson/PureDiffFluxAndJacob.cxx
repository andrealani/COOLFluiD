#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolumePoisson/PureDiffFluxAndJacob.hh"
#include "FiniteVolumePoisson/FiniteVolumePoisson.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "Poisson/PoissonDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<PureDiffFluxAndJacob,
                       CellCenterFVMData,
                       ComputeDiffusiveFlux,
                       FiniteVolumePoissonModule>
PureDiffFluxAndJacobProvider("PureDiffFluxAndJacob");

//////////////////////////////////////////////////////////////////////////////

PureDiffFluxAndJacob::PureDiffFluxAndJacob(const std::string& name) :
  PureDiffFlux(name)
{
}

//////////////////////////////////////////////////////////////////////////////

PureDiffFluxAndJacob::~PureDiffFluxAndJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void PureDiffFluxAndJacob::computeFlux(RealVector& result)
{
  if (!this->getMethodData().isPerturb()) {
    // reset the resulting flux to 0
    result = 0.0;
    
    GeometricEntity& geo = *getMethodData().getCurrentFace();
    SafePtr<DerivativeComputer> derivComputer = getMethodData().getDerivativeComputer();
    const bool isPerturb = this->getMethodData().isPerturb();

    // set the state values (pointers) corresponding to the vertices of the control volume
    derivComputer->computeControlVolume(_states, &geo);
    _nbCVStates = derivComputer->getNbVerticesInControlVolume(&geo);
    
    _varSet->setGradientVars(_states, _values, _nbCVStates);
    
    // compute control volume around the face and gradients
    derivComputer->computeGradients(&geo, _values, _gradients);
    
    // compute the average values
    derivComputer->computeAverageValues(&geo, _states, _avState);
    
    const CFreal faceArea = this->socket_faceAreas.getDataHandle()[geo.getID()];
    
    // set the flux
    result = _varSet->getFlux(_avState, _gradients, getMethodData().getUnitNormal());
    result *= faceArea;
    
    // jacobian computation
    SafePtr<vector<RealVector> > gradientsJacob = derivComputer->getGradientsJacob();
    RealVector& leftGradJacob  = (*gradientsJacob)[0];
    RealVector& rightGradJacob = (*gradientsJacob)[1];
    const CFreal radius = 0.;
    const RealVector& unitNormal = getMethodData().getUnitNormal();
    
    _varSet->computeFluxJacobian(_avState, leftGradJacob, unitNormal, radius, _lFluxJacobian);
    _lFluxJacobian *= faceArea; 
    
    if (!geo.getState(1)->isGhost()) {
      _varSet->computeFluxJacobian(_avState, rightGradJacob, unitNormal, radius, _rFluxJacobian);
      _rFluxJacobian *= faceArea;  
    }
  }
  else {
    PureDiffFlux::computeFlux(result);
  }
  
  CFLog(VERBOSE, "PureDiffFluxAndJacob::computeFlux() => result = " << result << "\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
