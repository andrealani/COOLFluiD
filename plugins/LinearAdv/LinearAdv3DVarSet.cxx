// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdv {


//////////////////////////////////////////////////////////////////////////////

void LinearAdv3DVarSet::computeFlux(const RealVector& pdata, const RealVector& normals)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];

  _fluxArray[0] = linearData[LinearAdvTerm::VX] *pdata[0]* nx +
    linearData[LinearAdvTerm::VY] *pdata[0]* ny + linearData[LinearAdvTerm::VZ] *pdata[0]* nz;
}


//////////////////////////////////////////////////////////////////////////////

void LinearAdv3DVarSet::computeStateFlux(const RealVector& pdata)
{
  // fx = vx*u, fy = vy*u
  _physFlux(0,XX) = pdata[LinearAdvTerm::VX];
  _physFlux(0,YY) = pdata[LinearAdvTerm::VY];
  _physFlux(0,ZZ) = pdata[LinearAdvTerm::VZ];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdv3DVarSet::getMaxEigenValue(const RealVector& pdata,
               const RealVector& normal){
  
  return pdata[LinearAdvTerm::VX]*normal[XX] + pdata[LinearAdvTerm::VY]*normal[YY] +pdata[LinearAdvTerm::VZ]*normal[ZZ];

}

//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdv3DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
            const RealVector& normal){

return fabs(pdata[LinearAdvTerm::VX]*normal[XX] + pdata[LinearAdvTerm::VY]*normal[YY] + pdata[LinearAdvTerm::VZ]*normal[ZZ]);

}

//////////////////////////////////////////////////////////////////////////////

void LinearAdv3DVarSet::computeEigenValues(const RealVector& pdata,
               const RealVector& normal, RealVector& result){

 result[0] = pdata[LinearAdvTerm::VX]*normal[XX] + pdata[LinearAdvTerm::VY]*normal[YY] + pdata[LinearAdvTerm::VZ]*normal[ZZ];

}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
