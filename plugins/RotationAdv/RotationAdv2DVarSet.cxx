// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdv2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {


//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DVarSet::computeFlux(const RealVector& pdata, const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  _fluxArray[0] = pdata[RotationAdvTerm::VX] * pdata[0] * nx + pdata[RotationAdvTerm::VY] * pdata[0] * ny;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DVarSet::computeStateFlux(const RealVector& pdata)
{
  // fx = vx*u, fy = vy*u
  _physFlux(0,XX) = pdata[RotationAdvTerm::VX]* pdata[0] ;
  _physFlux(0,YY) = pdata[RotationAdvTerm::VY]* pdata[0] ;
}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdv2DVarSet::getMaxEigenValue(const RealVector& pdata,
						 const RealVector& normal)
{
   return pdata[RotationAdvTerm::VX]*normal[XX] +
          pdata[RotationAdvTerm::VY]*normal[YY];
 
}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdv2DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
             const RealVector& normal)
{
   return fabs(pdata[RotationAdvTerm::VX]*normal[XX] +
          pdata[RotationAdvTerm::VY]*normal[YY]);

}

//////////////////////////////////////////////////////////////////////////////

void RotationAdv2DVarSet::computeEigenValues(const RealVector & pdata,
					     const RealVector& normal,
					     RealVector& result)
{
 
   result[0] = pdata[RotationAdvTerm::VX]*normal[XX] +
               pdata[RotationAdvTerm::VY]*normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
