// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "AdvectionDiffusionSys2DVarSet.hh"
#include "LinearAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


  namespace Physics {
namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusionSys2DVarSet::setup()
{
  AdvectionDiffusionSysVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& AdvectionDiffusionSys2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const RealVector& normal,
                                          const CFreal& radius)
{
  const RealVector& gradu = *gradients[0];
  RealVector& adData = getModel().getPhysicalData();
 adData[ADSysTerm::NU] = getModel().getDiffusionCoeff();

  _flux[0] = adData[ADSysTerm::NU]*(normal[XX]*gradu[XX] + normal[YY]*gradu[YY]);

 return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& AdvectionDiffusionSys2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
 setGradientState(state);
  const RealVector& gradu = *gradients[0];
  RealVector& adData = getModel().getPhysicalData();
  adData[ADSysTerm::NU] = getModel().getDiffusionCoeff();

  _fluxVec(0,XX) = adData[ADSysTerm::NU]*gradu[XX];
  _fluxVec(0,YY) = adData[ADSysTerm::NU]*gradu[YY];

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
