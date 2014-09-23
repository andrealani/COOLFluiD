// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdvSys/RotationAdvSys.hh"
#include "RotationAdvSys2DVarSet.hh"
#include "Environment/ObjectProvider.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

RotationAdvSys2DVarSet::RotationAdvSys2DVarSet(Common::SafePtr<BaseTerm> term) :
  RotationAdvSysVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvSys2DVarSet::~RotationAdvSys2DVarSet()
{
}

////////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DVarSet::setup()
{

   RotationAdvSysVarSet::setup();

//                  // set EquationSetData
   RotationAdvSys2DVarSet::getEqSetData().resize(1);
   RotationAdvSys2DVarSet::getEqSetData()[0].setup(0,0,4);
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DVarSet::computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues)
{
  throw Common::NotImplementedException (FromHere(),"");
}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdvSys2DVarSet::getMaxEigenValue(const RealVector& pdata, const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdvSys2DVarSet::getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal)
{
   throw Common::NotImplementedException (FromHere(),"");
   return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DVarSet::computeFlux(const RealVector& pdata, const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  _fluxArray[0] = pdata[RotationAdvSysTerm::C0X] * pdata[RotationAdvSysTerm::u0] * nx + pdata[RotationAdvSysTerm::C0Y] * pdata[RotationAdvSysTerm::u0] * ny;
  _fluxArray[1] = pdata[RotationAdvSysTerm::C1X] * pdata[RotationAdvSysTerm::u1] * nx + pdata[RotationAdvSysTerm::C1Y] * pdata[RotationAdvSysTerm::u1] * ny;
  _fluxArray[2] = pdata[RotationAdvSysTerm::C2X] * pdata[RotationAdvSysTerm::u2] * nx + pdata[RotationAdvSysTerm::C2Y] * pdata[RotationAdvSysTerm::u2] * ny;
  _fluxArray[3] = pdata[RotationAdvSysTerm::C3X] * pdata[RotationAdvSysTerm::u3] * nx + pdata[RotationAdvSysTerm::C3Y] * pdata[RotationAdvSysTerm::u3] * ny;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DVarSet::computeStateFlux(const RealVector& pdata)
{
  // fx = vx*u, fy = vy*u
  _physFlux(0,XX) = pdata[RotationAdvSysTerm::C0X]* pdata[RotationAdvSysTerm::u0] ;
  _physFlux(0,YY) = pdata[RotationAdvSysTerm::C0Y]* pdata[RotationAdvSysTerm::u0] ;
  
  _physFlux(1,XX) = pdata[RotationAdvSysTerm::C1X]* pdata[RotationAdvSysTerm::u1] ;
  _physFlux(1,YY) = pdata[RotationAdvSysTerm::C1Y]* pdata[RotationAdvSysTerm::u1] ;
  
  _physFlux(2,XX) = pdata[RotationAdvSysTerm::C2X]* pdata[RotationAdvSysTerm::u2] ;
  _physFlux(2,YY) = pdata[RotationAdvSysTerm::C2Y]* pdata[RotationAdvSysTerm::u2] ;
  
  _physFlux(3,XX) = pdata[RotationAdvSysTerm::C3X]* pdata[RotationAdvSysTerm::u3] ;
  _physFlux(3,YY) = pdata[RotationAdvSysTerm::C3Y]* pdata[RotationAdvSysTerm::u3] ;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
