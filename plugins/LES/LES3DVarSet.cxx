#include "LES/LES3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

extern "C" {

#include "les_interface.h"

} // extern "C"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace LES {

//////////////////////////////////////////////////////////////////////////////

void LES3DVarSet::setup()
{
  LESVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& LES3DVarSet::getFlux(const RealVector& state,
                                 const vector<RealVector*>& gradients,
                                 const RealVector& normal,
                                 const CFreal& radius)
{
  // set the gradient variables in the gradient state (not the gradients themselves!)
  setGradientState(state);

  // set pointer to current state
  m_currState = &state;

  // set pointer to current gradients (to enable passing the data to the library
  m_currGrad = &gradients;

  // set gradients and velocities
  const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const RealVector& gradW = *gradients[3];
  const RealVector& gradT = *gradients[4];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY] + gradW[ZZ]);
  const CFreal avU = _gradState[1];
  const CFreal avV = _gradState[2];
  const CFreal avW = _gradState[3];

  // get the vector to store the physical data in
  RealVector& nsData = getModel().getPhysicalData();

  // adimensional dynamic viscosity and thermal conductivity
  if (_useBackUpValues) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamic viscosity
    nsData[NSTerm::MU] = getDynViscosity(state, gradients)*getModel().getArtDiffCoeff();;

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] = getThermConductivity(state, nsData[NSTerm::MU]);

    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }

  // get adimensional eddy viscosity and thermal conductivity
  CFreal eddyDynVisc   = 0.0;
  les_compute_eddy_dynvisc(&eddyDynVisc  );
  CFreal eddyThermCond = 0.0;
  les_compute_thermal_cond(&eddyThermCond);

  // compute flux
  const CFreal coeffTauMu = getModel().getCoeffTau()*(nsData[NSTerm::MU] + eddyDynVisc);
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const CFreal tauXZ = coeffTauMu*(gradU[ZZ] + gradW[XX]);
  const CFreal tauYZ = coeffTauMu*(gradW[YY] + gradV[ZZ]);
  const CFreal tauZZ = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV);
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal qFlux = -getModel().getCoeffQ()*(nsData[NSTerm::LAMBDA] + eddyThermCond)*
                              (gradT[XX]*nx + gradT[YY]*ny + gradT[ZZ]*nz);

  _flux[1] = tauXX*nx + tauXY*ny + tauXZ*nz;
  _flux[2] = tauXY*nx + tauYY*ny + tauYZ*nz;
  _flux[3] = tauXZ*nx + tauYZ*ny + tauZZ*nz;
  _flux[4] = (tauXX*avU + tauXY*avV + tauXZ*avW)*nx +
             (tauXY*avU + tauYY*avV + tauYZ*avW)*ny +
             (tauXZ*avU + tauYZ*avV + tauZZ*avW)*nz - qFlux;

  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& LES3DVarSet::getFlux(const RealVector& state,
                                 const vector<RealVector*>& gradients,
                                 const CFreal& radius)
{
  throw Common::NotImplementedException (FromHere(),"RealMatrix& LES3DVarSet::getFlux() not implemented...");
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

CFreal LES3DVarSet::getHeatFlux(const RealVector& state,
                                const vector<RealVector*>& gradients,
                                const RealVector& normal)
{
  const CFreal mu = getDynViscosity(state, gradients);
  const CFreal lambda = getThermConductivity(state, mu);
  const RealVector& gradT = *gradients[4];

  // get adimensional eddy viscosity and thermal conductivity
  CFreal eddyDynVisc   = 0.0;
  les_compute_eddy_dynvisc(&eddyDynVisc  );
  CFreal eddyThermCond = 0.0;
  les_compute_thermal_cond(&eddyThermCond);

  return -getModel().getCoeffQ()*(lambda + eddyThermCond)*(gradT[XX]*normal[XX] + gradT[YY]*normal[YY] + gradT[ZZ]*normal[ZZ]);
}

//////////////////////////////////////////////////////////////////////////////

void LES3DVarSet::getAxiSourceTerm(const RealVector& physicalData,
                                   const RealVector& state,
                                   const vector<RealVector*>& gradients,
                                   const CFreal& radius,
                                   RealVector& source)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace LES

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
