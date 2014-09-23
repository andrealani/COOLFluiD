// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationDiffusion2DVarSet.hh"
#include "RotationAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

void RotationDiffusion2DVarSet::setup()
{
  RotationDiffusionVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& RotationDiffusion2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const RealVector& normal,
                                          const CFreal& radius)
{
  const RealVector& gradu = *gradients[0];
  RealVector& rdData = getModel().getPhysicalData();
 rdData[RDTerm::NU] = getModel().getDiffusionCoeff();

  _flux[0] = rdData[RDTerm::NU]*(normal[XX]*gradu[0] + normal[YY]*gradu[1]);

 return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& RotationDiffusion2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
 setGradientState(state);
  const RealVector& gradu = *gradients[0];
  RealVector& rdData = getModel().getPhysicalData();
  rdData[RDTerm::NU] = getModel().getDiffusionCoeff();

  _fluxVec(0,XX) = rdData[RDTerm::NU]*gradu[0];
  _fluxVec(0,YY) = rdData[RDTerm::NU]*gradu[0];

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
