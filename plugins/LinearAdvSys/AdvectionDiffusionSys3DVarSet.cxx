// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "AdvectionDiffusionSys3DVarSet.hh"
#include "LinearAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusionSys3DVarSet::setup()
{
  AdvectionDiffusionSysVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& AdvectionDiffusionSys3DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const RealVector& normal,
                                          const CFreal& radius)
{
  const RealVector& gradu = *gradients[0];
  RealVector& adData = getModel().getPhysicalData();
 adData[ADSysTerm::NU] = getModel().getDiffusionCoeff();

 for (CFuint i =0; i<4; i++)
   _flux[i] = adData[ADSysTerm::NU]*(normal[XX]*gradu[XX] + normal[YY]*gradu[YY] + normal[ZZ]*gradu[ZZ]);

 return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& AdvectionDiffusionSys3DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
 setGradientState(state);
  const RealVector& gradu = *gradients[0];
  RealVector& adData = getModel().getPhysicalData();
  adData[ADSysTerm::NU] = getModel().getDiffusionCoeff();
for (CFuint i =0; i<4; i++)
  {
  _fluxVec(i,XX) = adData[ADSysTerm::NU]*gradu[XX];
  _fluxVec(i,YY) = adData[ADSysTerm::NU]*gradu[YY];
  _fluxVec(i,ZZ) = adData[ADSysTerm::NU]*gradu[ZZ];
  }
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
