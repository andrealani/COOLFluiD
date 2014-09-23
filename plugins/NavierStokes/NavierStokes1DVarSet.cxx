#include "NavierStokes1DVarSet.hh"
#include "EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DVarSet::setup()
{
  NavierStokesVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes1DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const RealVector& normal,
                                          const CFreal& radius)
{
  setGradientState(state);

  const RealVector& gradU = *gradients[1];
  const RealVector& gradT = *gradients[2];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX]);
  const CFreal avU = _gradState[1];

  RealVector& nsData = getModel().getPhysicalData();

  // adimensional dynamical viscosity
  if (_useBackUpValues || _freezeDiffCoeff) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamical viscosity
    nsData[NSTerm::MU] = getDynViscosity(state, gradients)*getModel().getArtDiffCoeff();

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] = getThermConductivity
      (state, nsData[NSTerm::MU]);


    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }
  const CFreal coeffTauMu = getModel().getCoeffTau()*nsData[NSTerm::MU];
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal nx = normal[XX];
  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*(gradT[XX]*nx);

  _flux[1] = tauXX*nx;
  _flux[2] = tauXX*avU*nx - qFlux;

  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& NavierStokes1DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
  setGradientState(state);

  const RealVector& gradU = *gradients[1];
  const RealVector& gradT = *gradients[2];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX]);
  const CFreal avU = _gradState[1];

  RealVector& nsData = getModel().getPhysicalData();

  // adimensional dynamical viscosity
  if (_useBackUpValues || _freezeDiffCoeff) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamical viscosity
    nsData[NSTerm::MU] = getDynViscosity(state, gradients);

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] = getThermConductivity
      (state, nsData[NSTerm::MU]);


    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }
  const CFreal coeffTauMu = getModel().getCoeffTau()*nsData[NSTerm::MU];
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT;

  _fluxVec(1,XX) = tauXX;
  _fluxVec(2,XX) = tauXX*avU - qFlux[XX];

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes1DVarSet::getHeatFlux(const RealVector& state,
					 const vector<RealVector*>& gradients,
					 const RealVector& normal)
{
  const CFreal mu = getDynViscosity(state, gradients);
  const CFreal lambda = getThermConductivity(state, mu);
  const RealVector& gradT = *gradients[2];
  
  return (-getModel().getCoeffQ()*lambda*(gradT[XX]*normal[XX]));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
