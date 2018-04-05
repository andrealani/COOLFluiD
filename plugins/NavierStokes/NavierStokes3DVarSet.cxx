#include "NavierStokes3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DVarSet::setup()
{
  NavierStokesVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes3DVarSet::getFlux(const RealVector& state,
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
  _flux.slice(_uID, dim) = _tau*normal;
  _flux[_TID] = MathFunctions::innerProd(_flux.slice(_uID,dim), _gradState.slice(_uID,dim)) - qFlux;
  
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& NavierStokes3DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
  // AL: this needs to be adapted for 2D and 1/2
  cf_assert(PhysicalModelStack::getActive()->getDim() == 3);
  
  setGradientState(state);
  // dummy normal here
  computeTransportProperties(state, gradients, _normal);
  computeStressTensor(state, gradients, radius);
  
  RealVector& nsData = getModel().getPhysicalData();
  const RealVector& gradT = *gradients[_TID];
  const CFreal avU = _gradState[_uID];
  const CFreal avV = _gradState[_vID];
  const CFreal avW = _gradState[_wID];
  const CFreal tauXX = _tau(XX, XX);
  const CFreal tauXY = _tau(XX, YY);
  const CFreal tauYY = _tau(YY, YY);
  const CFreal tauXZ = _tau(XX, ZZ);
  const CFreal tauYZ = _tau(YY, ZZ);
  const CFreal tauZZ = _tau(ZZ, ZZ);
  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT;
  
  _fluxVec(_uID,XX) = tauXX;
  _fluxVec(_uID,YY) = tauXY;
  _fluxVec(_uID,ZZ) = tauXZ;
  
  _fluxVec(_vID,XX) = tauXY;
  _fluxVec(_vID,YY) = tauYY;
  _fluxVec(_vID,ZZ) = tauYZ;
  
  _fluxVec(_wID,XX) = tauXZ;
  _fluxVec(_wID,YY) = tauYZ;
  _fluxVec(_wID,ZZ) = tauZZ;
  
  _fluxVec(_TID,XX) = tauXX*avU + tauXY*avV + tauXZ*avW - qFlux[XX];
  _fluxVec(_TID,YY) = tauXY*avU + tauYY*avV + tauYZ*avW - qFlux[YY];
  _fluxVec(_TID,ZZ) = tauXZ*avU + tauYZ*avV + tauZZ*avW - qFlux[ZZ];
  
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
