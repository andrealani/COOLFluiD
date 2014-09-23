// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "LinearAdvSys/LinearAdvSys.hh"
#include "LinearAdvSys2DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

LinearAdvSys2DVarSet::LinearAdvSys2DVarSet(Common::SafePtr<BaseTerm> term) :
  LinearAdvSysVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

LinearAdvSys2DVarSet::~LinearAdvSys2DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DVarSet::setup()
{

/*
static int ImIn=0;
if(ImIn == 0)
	{
  cout<< "I`m in LinearAdvSys2DVarSet::setup() \n";
  cin.get();
  ImIn=1;
} 
*/
  LinearAdvSysVarSet::setup();

  // set EquationSetData
  LinearAdvSys2DVarSet::getEqSetData().resize(1);
  LinearAdvSys2DVarSet::getEqSetData()[0].setup(0,0,4);
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DVarSet::computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues)
{
  throw Common::NotImplementedException (FromHere(),"");
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdvSys2DVarSet::getMaxEigenValue(const RealVector& pdata, const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdvSys2DVarSet::getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DVarSet::computeFlux(const RealVector& pdata, const RealVector& normals)
{

 
   const RealVector& linearData = getModel()->getPhysicalData();
 
    const CFreal c0x = linearData[LinearAdvSysTerm::C0X];
    const CFreal c1x = linearData[LinearAdvSysTerm::C1X];
    const CFreal c2x = linearData[LinearAdvSysTerm::C2X];
    const CFreal c3x = linearData[LinearAdvSysTerm::C3X];

    const CFreal c0y = linearData[LinearAdvSysTerm::C0Y];
    const CFreal c1y = linearData[LinearAdvSysTerm::C1Y];
    const CFreal c2y = linearData[LinearAdvSysTerm::C2Y];
    const CFreal c3y = linearData[LinearAdvSysTerm::C3Y];


    const CFreal nx = normals[XX];
    const CFreal ny = normals[YY];
    const CFreal c0n = c0x*nx+c0y*ny;
    const CFreal c1n = c1x*nx+c1y*ny;
    const CFreal c2n = c2x*nx+c2y*ny;
    const CFreal c3n = c3x*nx+c3y*ny;


     _fluxArray[0] = c0n*pdata[LinearAdvSysTerm::u0];
     _fluxArray[1] = c1n*pdata[LinearAdvSysTerm::u1];
     _fluxArray[2] = c2n*pdata[LinearAdvSysTerm::u2];
     _fluxArray[3] = c3n*pdata[LinearAdvSysTerm::u3];



}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DVarSet::computeStateFlux(const RealVector& pdata)
{

const RealVector& linearData = getModel()->getPhysicalData();

  // fx = vx*u, fy = vy*u
  _physFlux(0,XX) = pdata[LinearAdvSysTerm::C0X]* pdata[LinearAdvSysTerm::u0] ;
  _physFlux(0,YY) = pdata[LinearAdvSysTerm::C0Y]* pdata[LinearAdvSysTerm::u0] ;
  _physFlux(1,XX) = pdata[LinearAdvSysTerm::C1X]* pdata[LinearAdvSysTerm::u1] ;
  _physFlux(1,YY) = pdata[LinearAdvSysTerm::C1Y]* pdata[LinearAdvSysTerm::u1] ;
  _physFlux(2,XX) = pdata[LinearAdvSysTerm::C2X]* pdata[LinearAdvSysTerm::u2] ;
  _physFlux(2,YY) = pdata[LinearAdvSysTerm::C2Y]* pdata[LinearAdvSysTerm::u2] ;
  _physFlux(3,XX) = pdata[LinearAdvSysTerm::C3X]* pdata[LinearAdvSysTerm::u3] ;
  _physFlux(3,YY) = pdata[LinearAdvSysTerm::C3Y]* pdata[LinearAdvSysTerm::u3] ;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdvSys

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
