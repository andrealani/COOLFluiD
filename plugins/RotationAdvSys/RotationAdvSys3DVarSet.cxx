// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdvSys3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {


RotationAdvSys3DVarSet::RotationAdvSys3DVarSet(Common::SafePtr<BaseTerm> term) :
  RotationAdvSysVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvSys3DVarSet::~RotationAdvSys3DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////////

void RotationAdvSys3DVarSet::setup()
{

   RotationAdvSysVarSet::setup();
   // set EquationSetData
   RotationAdvSys3DVarSet::getEqSetData().resize(1);
   RotationAdvSys3DVarSet::getEqSetData()[0].setup(0,0,4);
}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys3DVarSet::computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues)
{
 throw Common::NotImplementedException (FromHere(),"");
}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdvSys3DVarSet::getMaxEigenValue(const RealVector& pdata, const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdvSys3DVarSet::getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"");
  return 0.;
}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys3DVarSet::computeFlux(const RealVector& pdata, const RealVector& normals)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];

  _fluxArray[0] = linearData[RotationAdvSysTerm::C0X] *pdata[0]* nx +
                  linearData[RotationAdvSysTerm::C0Y] *pdata[0]* ny +
                  linearData[RotationAdvSysTerm::C0Z] *pdata[0]* nz;

 _fluxArray[1] = linearData[RotationAdvSysTerm::C1X] *pdata[1]* nx +    
                 linearData[RotationAdvSysTerm::C1Y] *pdata[1]* ny + 
                 linearData[RotationAdvSysTerm::C1Z] *pdata[1]* nz;

 _fluxArray[2] = linearData[RotationAdvSysTerm::C2X] *pdata[2]* nx +    
                 linearData[RotationAdvSysTerm::C2Y] *pdata[2]* ny + 
                 linearData[RotationAdvSysTerm::C2Z] *pdata[2]* nz;

 _fluxArray[3] = linearData[RotationAdvSysTerm::C3X] *pdata[3]* nx +    
                 linearData[RotationAdvSysTerm::C3Y] *pdata[3]* ny + 
                 linearData[RotationAdvSysTerm::C3Z] *pdata[3]* nz;


}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys3DVarSet::computeStateFlux(const RealVector& pdata)
{

 const RealVector& linearData = getModel()->getPhysicalData();

// fx = vx*u, fy = vy*u
  _physFlux(0,XX) = pdata[RotationAdvSysTerm::C0X]*pdata[0];
  _physFlux(0,YY) = pdata[RotationAdvSysTerm::C0Y]*pdata[0];
  _physFlux(0,ZZ) = pdata[RotationAdvSysTerm::C0Z]*pdata[0];

  _physFlux(1,XX) = pdata[RotationAdvSysTerm::C1X]*pdata[1];
  _physFlux(1,YY) = pdata[RotationAdvSysTerm::C1Y]*pdata[1];
  _physFlux(1,ZZ) = pdata[RotationAdvSysTerm::C1Z]*pdata[1];

  _physFlux(2,XX) = pdata[RotationAdvSysTerm::C2X]*pdata[2];
  _physFlux(2,YY) = pdata[RotationAdvSysTerm::C2Y]*pdata[2];
  _physFlux(2,ZZ) = pdata[RotationAdvSysTerm::C2Z]*pdata[2];

  _physFlux(3,XX) = pdata[RotationAdvSysTerm::C3X]*pdata[3];
  _physFlux(3,YY) = pdata[RotationAdvSysTerm::C3Y]*pdata[3];
  _physFlux(3,ZZ) = pdata[RotationAdvSysTerm::C3Z]*pdata[3];


}



//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
