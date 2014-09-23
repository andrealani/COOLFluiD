// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RotationAdvSys/RotationAdvSys.hh"
#include "RotationAdvSys/RotationAdvSys2DPrim.hh"
#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

//////////////////////////////////////////// //////////////////////////

Environment::ObjectProvider<RotationAdvSys2DPrim, ConvectiveVarSet, RotationAdvSysModule, 1>
rotationAdvSys2DPrimProvider("RotationAdvSys2DPrim");

//////////////////////////////////////////////////////////////////////////////

RotationAdvSys2DPrim::RotationAdvSys2DPrim(Common::SafePtr<BaseTerm> term) :
  RotationAdvSys2DVarSet(term),
  _rightEv(4,4),
  _leftEv(4,4)
{
  vector<std::string> names(4);
  names[0] = "u0";
  names[1] = "u1";
  names[2] = "u2";
  names[3] = "u3";

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

RotationAdvSys2DPrim::~RotationAdvSys2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::setup()
{
  RotationAdvSys2DVarSet::setup();

    setConstJacob();
}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::setConstJacob()
{

 vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[XX](0,1) = 0.0;
  (*jacobians)[XX](0,2) = 0.0;
  (*jacobians)[XX](0,3) = 0.0;

  (*jacobians)[XX](1,0) = 0.0;
  (*jacobians)[XX](1,2) = 0.0;
  (*jacobians)[XX](1,3) = 0.0;

  (*jacobians)[XX](2,0) = 0.0;
  (*jacobians)[XX](2,1) = 0.0;
  (*jacobians)[XX](2,3) = 0.0;

  (*jacobians)[XX](3,0) = 0.0;
  (*jacobians)[XX](3,1) = 0.0;
  (*jacobians)[XX](3,2) = 0.0;



  (*jacobians)[YY](0,1) = 0.0;
  (*jacobians)[YY](0,2) = 0.0;
  (*jacobians)[YY](0,3) = 0.0;

  (*jacobians)[YY](1,0) = 0.0;
  (*jacobians)[YY](1,2) = 0.0;
  (*jacobians)[YY](1,3) = 0.0;

  (*jacobians)[YY](2,0) = 0.0;
  (*jacobians)[YY](2,1) = 0.0;
  (*jacobians)[YY](2,3) = 0.0;

  (*jacobians)[YY](3,0) = 0.0;
  (*jacobians)[YY](3,1) = 0.0;
  (*jacobians)[YY](3,2) = 0.0;


}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::computeProjectedJacobian(const RealVector& normal,
                                           RealMatrix& jacob)

{
  const RealVector& linearData = getModel()->getPhysicalData();

	
  jacob(0,0) = linearData[RotationAdvSysTerm::C0X]*normal[XX]+linearData[RotationAdvSysTerm::C0Y]*normal[YY];
  jacob(0,1) = 0.0;
  jacob(0,2) = 0.0;
  jacob(0,3) = 0.;	
	
  jacob(1,0) = 0.;
  jacob(1,1) = linearData[RotationAdvSysTerm::C1X]*normal[XX]+linearData[RotationAdvSysTerm::C1Y]*normal[YY];
  jacob(1,2) = 0.;
  jacob(1,3) = 0.;

  jacob(2,0) = 0.;
  jacob(2,1) = 0.;
  jacob(2,2) = linearData[RotationAdvSysTerm::C2X]*normal[XX]+linearData[RotationAdvSysTerm::C2Y]*normal[YY];
  jacob(2,3) = 0.;

  jacob(3,0) = 0.;
  jacob(3,1) = 0.;
  jacob(3,2) = 0.;
  jacob(3,3) = linearData[RotationAdvSysTerm::C3X]*normal[XX]+linearData[RotationAdvSysTerm::C3Y]*normal[YY];
	
}

//////////////////////////////////////////////////////////////////////////////
//
void RotationAdvSys2DPrim::computeJacobians()
{
  const RealVector& lineardata = getModel()->getPhysicalData();
  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();
 
        (*jacobians)[XX](0,0) = lineardata[RotationAdvSysTerm::C0X];
	(*jacobians)[XX](1,1) = lineardata[RotationAdvSysTerm::C1X];
        (*jacobians)[XX](2,2) = lineardata[RotationAdvSysTerm::C2X];
        (*jacobians)[XX](2,2) = lineardata[RotationAdvSysTerm::C3X];

        (*jacobians)[YY](0,0) = lineardata[RotationAdvSysTerm::C0Y];
        (*jacobians)[YY](1,1) = lineardata[RotationAdvSysTerm::C1Y];
        (*jacobians)[YY](2,2) = lineardata[RotationAdvSysTerm::C2Y];
        (*jacobians)[YY](2,2) = lineardata[RotationAdvSysTerm::C3Y];
  
}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
            RealMatrix& leftEv,
            RealVector& eValues,
            const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal c0x = linearData[RotationAdvSysTerm::C0X];
  const CFreal c1x = linearData[RotationAdvSysTerm::C1X];
  const CFreal c2x = linearData[RotationAdvSysTerm::C2X];
  const CFreal c3x = linearData[RotationAdvSysTerm::C3X];

  const CFreal c0y = linearData[RotationAdvSysTerm::C0Y];
  const CFreal c1y = linearData[RotationAdvSysTerm::C1Y];
  const CFreal c2y = linearData[RotationAdvSysTerm::C2Y];
  const CFreal c3y = linearData[RotationAdvSysTerm::C3Y];

  const CFreal Vn0 = c0x*nx + c0y*ny;
  const CFreal Vn1 = c1x*nx + c1y*ny;
  const CFreal Vn2 = c2x*nx + c2y*ny;
  const CFreal Vn3 = c3x*nx + c3y*ny;


 rightEv(0,0) = 1.0;  rightEv(1,0) = 0.0;  rightEv(2,0) = 0.0;  rightEv(3,0) = 0.0;
 rightEv(0,1) = 0.0;  rightEv(1,1) = 1.0;  rightEv(2,1) = 0.0;  rightEv(3,1) = 0.0;
 rightEv(0,2) = 0.0;  rightEv(1,2) = 0.0;  rightEv(2,2) = 1.0;  rightEv(3,2) = 0.0;
 rightEv(0,3) = 0.0;  rightEv(1,3) = 0.0;  rightEv(2,3) = 0.0;  rightEv(3,3) = 1.0;


  leftEv(0,0) = 1.0; leftEv(0,1) = 0.0; leftEv(0,2) = 0.0; leftEv(0,3) = 0.0; 
  leftEv(1,0) = 0.0; leftEv(1,1) = 1.0; leftEv(1,2) = 0.0; leftEv(1,3) = 0.0;
  leftEv(2,0) = 0.0; leftEv(2,1) = 0.0; leftEv(2,2) = 1.0; leftEv(2,3) = 0.0;
  leftEv(3,0) = 0.0; leftEv(3,1) = 0.0; leftEv(3,2) = 0.0; leftEv(3,3) = 1.0;

  eValues[0] = Vn0;
 eValues[1] = Vn1;
 eValues[2] = Vn2;
 eValues[3] = Vn3;
}



//////////////////////////////////////////////////////////////////////////////
//
void RotationAdvSys2DPrim::setEigenVect1(RealVector& r1,
                                       State& state,
                                       const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"RotationAdvSys2DPrim::setEigenVect1");
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::setEigenVect2(RealVector& r2,
                                       State& state,
                                       const RealVector& normal)
{
 throw Common::NotImplementedException (FromHere(),"RotationAdvSys2DPrim::setEigenVect2");
}
//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::setEigenVect3(RealVector& r3,
                                       State& state,
                                       const RealVector& normal)
{
 throw Common::NotImplementedException (FromHere(),"RotationAdvSys2DPrim::setEigenVect3");
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::setEigenVect4(RealVector& r4,
                                       State& state,
                                       const RealVector& normal)
{
 throw Common::NotImplementedException (FromHere(),"RotationAdvSys2DPrim::setEigenVect4");
}


//////////////////////////////////////////////////////////////////////////////

CFuint RotationAdvSys2DPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
 }

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  
  const CFreal c0x = linearData[RotationAdvSysTerm::C0X];
  const CFreal c1x = linearData[RotationAdvSysTerm::C1X];
  const CFreal c2x = linearData[RotationAdvSysTerm::C2X];
  const CFreal c3x = linearData[RotationAdvSysTerm::C3X];

  const CFreal c0y = linearData[RotationAdvSysTerm::C0Y];
  const CFreal c1y = linearData[RotationAdvSysTerm::C1Y];
  const CFreal c2y = linearData[RotationAdvSysTerm::C2Y];
  const CFreal c3y = linearData[RotationAdvSysTerm::C3Y];


  const CFreal Vn0 = c0x*nx + c0y*ny;
  const CFreal Vn1 = c1x*nx + c1y*ny;
  const CFreal Vn2 = c2x*nx + c2y*ny;
  const CFreal Vn3 = c3x*nx + c3y*ny;


 _rightEv(0,0) = 1.0;  _rightEv(1,0) = 0.0;  _rightEv(2,0) = 0.0;  _rightEv(3,0) = 0.0;
 _rightEv(0,1) = 0.0;  _rightEv(1,1) = 1.0;  _rightEv(2,1) = 0.0;  _rightEv(3,1) = 0.0;
 _rightEv(0,2) = 0.0;  _rightEv(1,2) = 0.0;  _rightEv(2,2) = 1.0;  _rightEv(3,2) = 0.0;
 _rightEv(0,3) = 0.0;  _rightEv(1,3) = 0.0;  _rightEv(2,3) = 0.0;  _rightEv(3,3) = 1.0;


  _leftEv(0,0) = 1.0; _leftEv(0,1) = 0.0; _leftEv(0,2) = 0.0; _leftEv(0,3) = 0.0;
  _leftEv(1,0) = 0.0; _leftEv(1,1) = 1.0; _leftEv(1,2) = 0.0; _leftEv(1,3) = 0.0;
  _leftEv(2,0) = 0.0; _leftEv(2,1) = 0.0; _leftEv(2,2) = 1.0; _leftEv(2,3) = 0.0;
  _leftEv(3,0) = 0.0; _leftEv(3,1) = 0.0; _leftEv(3,2) = 0.0; _leftEv(3,3) = 1.0;

  eValues[0] = Vn0;
 eValues[1] = Vn1;
 eValues[2] = Vn2;
 eValues[3] = Vn3;

  for (CFuint iEq = 0; iEq < eValues.size(); ++iEq)
  {
      _eValuesP[iEq] = max(0.,eValues[iEq]);
      _eValuesM[iEq] = min(0.,eValues[iEq]);
  }

  // compute jacobian + and -
     jacobPlus = _rightEv*(_eValuesP*_leftEv);
     jacobMin  = _rightEv*(_eValuesM*_leftEv);
  
}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::computePhysicalData(const State& state,
					    RealVector& data)
{
  RealVector& node = state.getCoordinates();

  const CFreal OX0 = getModel()->getOX0();
  const CFreal OY0 = getModel()->getOY0();

  const CFreal OX1 = getModel()->getOX1();
  const CFreal OY1 = getModel()->getOY1();

  const CFreal OX2 = getModel()->getOX2();
  const CFreal OY2 = getModel()->getOY2();

  const CFreal OX3 = getModel()->getOX3();
  const CFreal OY3 = getModel()->getOY3();


  data[RotationAdvSysTerm::u0] = state[0]; 
  data[RotationAdvSysTerm::u1] = state[1];
  data[RotationAdvSysTerm::u2] = state[2];
  data[RotationAdvSysTerm::u3] = state[3];


  data[RotationAdvSysTerm::C0X] = - (node[YY] - OY0);
  data[RotationAdvSysTerm::C0Y] =   (node[XX] - OX0);

  data[RotationAdvSysTerm::C1X] = - (node[YY] - OY1);
  data[RotationAdvSysTerm::C1Y] =   (node[XX] - OX1);

  data[RotationAdvSysTerm::C2X] = - (node[YY] - OY2);
  data[RotationAdvSysTerm::C2Y] =   (node[XX] - OX2);

  data[RotationAdvSysTerm::C3X] = - (node[YY] - OY3);
  data[RotationAdvSysTerm::C3Y] =   (node[XX] - OX3);


  if (getModel()->isClockwise()) data *= -1.0;
}



//////////////////////////////////////////////////////////////////////////////

CFreal RotationAdvSys2DPrim::getSpeed(const State& state) const
{
  throw Common::NotImplementedException (FromHere(),"RotationAdvSys2DPrim::getSpeed");
  return 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::setDimensionalValues(const State& state, RealVector& result)
{
   result = state;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::setAdimensionalValues(const State& state, RealVector& result)
{
    result = state;
}

//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::computePerturbedStatesData(const vector<State*>& states,
                                                    const CFuint nbStatesInVec,
                                                    const CFuint iVar)
{
   throw Common::NotImplementedException (FromHere(),"LinearAdvSys2DPrim::computePerturbedStatesData");
}


//////////////////////////////////////////////////////////////////////////////

void RotationAdvSys2DPrim::computeStateFromPhysicalData(const RealVector& data,
                                              State& state)
{
  state[0] = data[RotationAdvSysTerm::u0];
  state[1] = data[RotationAdvSysTerm::u1];
  state[2] = data[RotationAdvSysTerm::u2];
  state[3] = data[RotationAdvSysTerm::u3];	


}
                                                


//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
