// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdv2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdv {


//////////////////////////////////////////////////////////////////////////////

void LinearAdv2DVarSet::computeFlux(const RealVector& pdata, const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];

  cf_assert(pdata.size() > LinearAdvTerm::VY);
  
  _fluxArray[0] = pdata[LinearAdvTerm::VX] * pdata[0] * nx +
    pdata[LinearAdvTerm::VY] * pdata[0] * ny;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdv2DVarSet::computeStateFlux(const RealVector& pdata)
{
  // fx = vx*u, fy = vy*u
  _physFlux(0,XX) = pdata[LinearAdvTerm::VX]* pdata[0] ;
  _physFlux(0,YY) = pdata[LinearAdvTerm::VY]* pdata[0] ;
}


//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdv2DVarSet::getMaxEigenValue(const RealVector& pdata,
               const RealVector& normal){
  
  return pdata[LinearAdvTerm::VX]*normal[XX] + pdata[LinearAdvTerm::VY]*normal[YY];

}

//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdv2DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
            const RealVector& normal){

return fabs(pdata[LinearAdvTerm::VX]*normal[XX] + pdata[LinearAdvTerm::VY]*normal[YY]);

}

//////////////////////////////////////////////////////////////////////////////

void LinearAdv2DVarSet::computeEigenValues(const RealVector& pdata,
               const RealVector& normal, RealVector& result){

 result[0] = pdata[LinearAdvTerm::VX]*normal[XX] + pdata[LinearAdvTerm::VY]*normal[YY];

}


//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
