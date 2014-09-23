#include "NavierStokes2DVarSet.hh"
#include "EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DVarSet::setup()
{
  NavierStokesVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& NavierStokes2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
  setGradientState(state);
  // dummy normal here
  computeTransportProperties(state, gradients, _normal);
  computeStressTensor(state, gradients, radius);

  RealVector& nsData = getModel().getPhysicalData();
  const RealVector& gradT = *gradients[_TID];
  const CFreal avU = _gradState[_uID];
  const CFreal avV = _gradState[_vID];
  const CFreal tauXX = _tau(XX, XX);
  const CFreal tauXY = _tau(XX, YY);
  const CFreal tauYY = _tau(YY, YY);
  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT;
  
  _fluxVec(_uID,XX) = tauXX;
  _fluxVec(_uID,YY) = tauXY;
  
  _fluxVec(_vID,XX) = tauXY;
  _fluxVec(_vID,YY) = tauYY;
  
  _fluxVec(_TID,XX) = tauXX*avU + tauXY*avV - qFlux[XX];
  _fluxVec(_TID,YY) = tauXY*avU + tauYY*avV - qFlux[YY];
  
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DVarSet::getAxiSourceTerm(const RealVector& physicalData,
					    const RealVector& state,
					    const vector<RealVector*>& gradients,
					    const CFreal& radius,
					    RealVector& source)
{
  const CFreal rho = physicalData[EulerTerm::RHO];
  const CFreal u = physicalData[EulerTerm::VX];
  const CFreal v = physicalData[EulerTerm::VY];
  const CFreal rhov = rho*v;
  
  const CFreal vOverRadius = v/radius;
  const RealVector& gradU = *gradients[_uID];
  const RealVector& gradV = *gradients[_vID];
  const RealVector& gradT = *gradients[_TID];
  
  const CFreal twoThird = 2./3.;
  const CFreal divV = gradU[XX] + gradV[YY];
  const CFreal mu = getDynViscosity(state, gradients);
  const CFreal lambda = getThermConductivity(state, mu);
  const CFreal coeffTauMu = getModel().getCoeffTau()*mu;
  const CFreal tauRX = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauRR = coeffTauMu*(2.*gradV[YY] - twoThird*(divV + vOverRadius));
  const CFreal tauTT = -coeffTauMu*twoThird*(divV - 2.*vOverRadius);
  const CFreal qr = -getModel().getCoeffQ()*lambda*gradT[YY];
  
  source[_uID-1] = -rhov;
  source[_uID] = -rhov*u + tauRX;
  source[_vID] = -rhov*v + tauRR - tauTT;
  source[_TID] = -rhov*physicalData[EulerTerm::H] - qr + tauRX*u + tauRR*v;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes2DVarSet::getFlux(const RealVector& state,
					  const vector<RealVector*>& gradients,
					  const RealVector& normal,
					  const CFreal& radius)
{  
  setGradientState(state);
  computeTransportProperties(state, gradients, normal);
  computeStressTensor(state, gradients, radius);
  
  const RealVector& nsData = getModel().getPhysicalData();
  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*
    MathFunctions::innerProd(*gradients[_TID], normal);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _flux.slice(_uID,dim) = _tau*normal;
  // _flux[_TID] = MathFunctions::innerProd(_flux.slice(_uID,dim), _gradState.slice(_uID,dim)) - qFlux;
  
  // AL: old implementation, slower 
  _flux[_TID] = (_tau(XX,XX)*_gradState[_uID] + _tau(XX,YY)*_gradState[_vID])*normal[XX] + 
    (_tau(XX,YY)*_gradState[_uID] + _tau(YY,YY)*_gradState[_vID])*normal[YY] - qFlux;
    
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
