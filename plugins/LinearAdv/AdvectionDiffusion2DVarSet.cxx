// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "AdvectionDiffusion2DVarSet.hh"
#include "LinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

void AdvectionDiffusion2DVarSet::setup()
{
  AdvectionDiffusionVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& AdvectionDiffusion2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const RealVector& normal,
                                          const CFreal& radius)
{
  const RealVector& gradu = *gradients[0];
  RealVector& adData = getModel().getPhysicalData();
 adData[ADTerm::NU] = getModel().getDiffusionCoeff();

  _flux[0] = adData[ADTerm::NU]*(normal[XX]*gradu[XX] + normal[YY]*gradu[YY]);

 return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& AdvectionDiffusion2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
 setGradientState(state);
  const RealVector& gradu = *gradients[0];
  RealVector& adData = getModel().getPhysicalData();
  adData[ADTerm::NU] = getModel().getDiffusionCoeff();

  _fluxVec(0,XX) = adData[ADTerm::NU]*gradu[XX];
  _fluxVec(0,YY) = adData[ADTerm::NU]*gradu[YY];

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
