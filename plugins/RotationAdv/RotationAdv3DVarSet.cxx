// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////
void RotationAdv3DVarSet::computeFlux(const RealVector& pdata, const RealVector& normals)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];

  _fluxArray[0] = linearData[RotationAdvTerm::VX] *pdata[0]* nx +
    linearData[RotationAdvTerm::VY] *pdata[0]* ny + linearData[RotationAdvTerm::VZ] *pdata[0]* nz;
}


//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DVarSet::computeStateFlux(const RealVector& pdata)
{
  // fx = vx*u, fy = vy*u
  _physFlux(0,XX) = pdata[RotationAdvTerm::VX]*pdata[0];
  _physFlux(0,YY) = pdata[RotationAdvTerm::VY]*pdata[0];
  _physFlux(0,ZZ) = pdata[RotationAdvTerm::VZ]*pdata[0];
}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdv3DVarSet::getMaxEigenValue(const RealVector& pdata,
                                           const RealVector& normal){

  return pdata[RotationAdvTerm::VX]*normal[XX] + pdata[RotationAdvTerm::VY]*normal[YY] +pdata[RotationAdvTerm::VZ]*normal[ZZ];

}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdv3DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
                                              const RealVector& normal){

   return fabs(pdata[RotationAdvTerm::VX]*normal[XX] + pdata[RotationAdvTerm::VY]*normal[YY] + pdata[RotationAdvTerm::VZ]*normal[ZZ]);

}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv3DVarSet::computeEigenValues(const RealVector& pdata,
                                           const RealVector& normal, RealVector& result){

     result[0] = pdata[RotationAdvTerm::VX]*normal[XX] + pdata[RotationAdvTerm::VY]*normal[YY] + pdata[RotationAdvTerm::VZ]*normal[ZZ];

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
