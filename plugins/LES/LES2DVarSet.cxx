#include "LES/LES2DVarSet.hh"
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

void LES2DVarSet::setup()
{
  LESVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& LES2DVarSet::getFlux(const RealVector& state,
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
  const RealVector& gradT = *gradients[3];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY]);
  const CFreal avU = _gradState[1];
  const CFreal avV = _gradState[2];

  // get the vector to store the physical data in
  RealVector& nsData = getModel().getPhysicalData();

  // adimensional dynamic viscosity and thermal conductivity
  if (_useBackUpValues || _freezeDiffCoeff) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamic viscosity
    nsData[NSTerm::MU] = getDynViscosity(state, gradients)*getModel().getArtDiffCoeff();;

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] = getThermConductivity
      (state, nsData[NSTerm::MU]);


    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }

  // get adimensional eddy viscosity and thermal conductivity
  CFreal eddyDynVisc   = 0.0;
  les_compute_eddy_dynvisc(&eddyDynVisc  );
  CFreal eddyThermCond = 0.0;
  les_compute_thermal_cond(&eddyThermCond); // positive

  // compute flux
  const CFreal coeffTauMu = getModel().getCoeffTau()*(nsData[NSTerm::MU] + eddyDynVisc);
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal qFlux = -getModel().getCoeffQ()*(nsData[NSTerm::LAMBDA] + eddyThermCond)*
                        (gradT[XX]*nx + gradT[YY]*ny); // negative

  _flux[1] = tauXX*nx + tauXY*ny;
  _flux[2] = tauXY*nx + tauYY*ny;
  _flux[3] = (tauXX*avU + tauXY*avV)*nx + (tauXY*avU + tauYY*avV)*ny - qFlux;

  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& LES2DVarSet::getFlux(const RealVector& state,
                                 const vector<RealVector*>& gradients,
                                 const CFreal& radius)
{
/*  setGradientState(state);

  CFreal vOverRadius = 0.0;
  if (radius > MathTools::MathConsts::CFrealEps()) {
    // if the face is a boundary face, the radius could be 0
    // check against eps instead of 0. for safety
    vOverRadius = _gradState[2]/radius;
  }

  const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const RealVector& gradT = *gradients[3];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY] + vOverRadius);
  const CFreal avU = _gradState[1];
  const CFreal avV = _gradState[2];

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
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT;

  _fluxVec(1,XX) = tauXX;
  _fluxVec(1,YY) = tauXY;

  _fluxVec(2,XX) = tauXY;
  _fluxVec(2,YY) = tauYY;

  _fluxVec(3,XX) = tauXX*avU + tauXY*avV - qFlux[XX];
  _fluxVec(3,YY) = tauXY*avU + tauYY*avV - qFlux[YY];
*/
  throw Common::NotImplementedException (FromHere(),"RealMatrix& LES2DVarSet::getFlux() not implemented...");
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

CFreal LES2DVarSet::getHeatFlux(const RealVector& state,
                                const vector<RealVector*>& gradients,
                                const RealVector& normal)
{
  const CFreal mu = getDynViscosity(state, gradients);
  const CFreal lambda = getThermConductivity(state, mu);
  const RealVector& gradT = *gradients[3];

  // get adimensional eddy viscosity and thermal conductivity
  CFreal eddyDynVisc   = 0.0;
  les_compute_eddy_dynvisc(&eddyDynVisc  );
  CFreal eddyThermCond = 0.0;
  les_compute_thermal_cond(&eddyThermCond);

  return -getModel().getCoeffQ()*(lambda + eddyThermCond)*(gradT[XX]*normal[XX] + gradT[YY]*normal[YY]);
}

//////////////////////////////////////////////////////////////////////////////

void LES2DVarSet::getAxiSourceTerm(const RealVector& physicalData,
                                   const RealVector& state,
                                   const vector<RealVector*>& gradients,
                                   const CFreal& radius,
                                   RealVector& source)
{
    throw Common::NotImplementedException (FromHere(),"LES2DVarSet::getAxiSourceTerm()");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace LES

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
