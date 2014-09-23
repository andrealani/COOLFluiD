// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "LinearAdvSys/LinearAdvSys.hh"
#include "LinearAdvSys/LinearAdvSys2DPrim.hh"
#include "MathTools/MatrixInverter.hh"
#include "LinearAdvSys/LinearAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinearAdvSys2DPrim, ConvectiveVarSet, LinearAdvSysModule, 1>
linearAdvSys2DPrimProvider("LinearAdvSys2DPrim");

//////////////////////////////////////////////////////////////////////////////

LinearAdvSys2DPrim::LinearAdvSys2DPrim(Common::SafePtr<BaseTerm> term) :
  LinearAdvSys2DVarSet(term),
  _rightEv(4,4),
  _leftEv(4,4)
{
  ///@todo NV:Need to check if the conservative variable are that or something else

  vector<std::string> names(4);
///@Tom
  names[0] = "u0";
  names[1] = "u1";
  names[2] = "u2";
  names[3] = "u3";

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

LinearAdvSys2DPrim::~LinearAdvSys2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setup()
{
  LinearAdvSys2DVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setConstJacob()
{
 
 vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  Common::SafePtr<LinearAdvSysTerm> term = getModel();

  const CFreal c0x = term->getc0x();
  const CFreal c1x = term->getc1x();
  const CFreal c2x = term->getc2x();
  const CFreal c3x = term->getc3x();

  const CFreal c0y = term->getc0y();
  const CFreal c1y = term->getc1y();
  const CFreal c2y = term->getc2y();
  const CFreal c3y = term->getc3y();

  (*jacobians)[XX](0,0) = c0x;
  (*jacobians)[XX](0,1) = 0.0;
  (*jacobians)[XX](0,2) = 0.0;
  (*jacobians)[XX](0,3) = 0.0;

  (*jacobians)[XX](1,0) = 0.0;
  (*jacobians)[XX](1,1) = c1x;
  (*jacobians)[XX](1,2) = 0.0;
  (*jacobians)[XX](1,3) = 0.0;

  (*jacobians)[XX](2,0) = 0.0;
  (*jacobians)[XX](2,1) = 0.0;
  (*jacobians)[XX](2,2) = c2x;
  (*jacobians)[XX](2,3) = 0.0;

  (*jacobians)[XX](3,0) = 0.0;
  (*jacobians)[XX](3,1) = 0.0;
  (*jacobians)[XX](3,2) = 0.0;
  (*jacobians)[XX](3,3) = c3x;



  (*jacobians)[YY](0,0) = c0y;
  (*jacobians)[YY](0,1) = 0.0;
  (*jacobians)[YY](0,2) = 0.0;
  (*jacobians)[YY](0,3) = 0.0;

  (*jacobians)[YY](1,0) = 0.0;
  (*jacobians)[YY](1,1) = c1y;
  (*jacobians)[YY](1,2) = 0.0;
  (*jacobians)[YY](1,3) = 0.0;

  (*jacobians)[YY](2,0) = 0.0;
  (*jacobians)[YY](2,1) = 0.0;
  (*jacobians)[YY](2,2) = c2y;
  (*jacobians)[YY](2,3) = 0.0;

  (*jacobians)[YY](3,0) = 0.0;
  (*jacobians)[YY](3,1) = 0.0;
  (*jacobians)[YY](3,2) = 0.0;
  (*jacobians)[YY](3,3) = c3y;

}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob)
{
  /// @note TK: Because of the Jacobian is given by main flow it is constnt
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::computeJacobians()
{
  /// @note TK: Because of the Jacobian is given by main flow it is constnt
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
          RealMatrix& leftEv,
          RealVector& eValues,
          const RealVector& normal)
{
  Common::SafePtr<LinearAdvSysTerm> term = getModel();
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal c0x = term->getc0x();
  const CFreal c1x = term->getc1x();
  const CFreal c2x = term->getc2x();
  const CFreal c3x = term->getc3x();

  const CFreal c0y = term->getc0y();
  const CFreal c1y = term->getc1y();
  const CFreal c2y = term->getc2y();
  const CFreal c3y = term->getc3y();

  const CFreal Vn0 = c0x*nx + c0y*ny;
  const CFreal Vn1 = c1x*nx + c1y*ny;
  const CFreal Vn2 = c2x*nx + c2y*ny;
  const CFreal Vn3 = c3x*nx + c3y*ny;

  rightEv(0,0) = 1.0;
  rightEv(0,1) = 0.0;
  rightEv(0,2) = 0.0;
  rightEv(0,3) = 0.0;

  rightEv(1,0) = 0.0;
  rightEv(1,1) = 1.0;
  rightEv(1,2) = 0.0;
  rightEv(1,3) = 0.0;

  rightEv(2,0) = 0.0;
  rightEv(2,1) = 0.0;
  rightEv(2,2) = 1.0;
  rightEv(2,3) = 0.0;

  rightEv(3,0) = 0.0;
  rightEv(3,1) = 0.0;
  rightEv(3,2) = 0.0;
  rightEv(3,3) = 1.0;



  leftEv(0,0) = 1.0;
  leftEv(0,1) = 0.0;
  leftEv(0,2) = 0.0;
  leftEv(0,3) = 0.0;

  leftEv(1,0) = 0.0;
  leftEv(1,1) = 1.0;
  leftEv(1,2) = 0.0;
  leftEv(1,3) = 0.0;

  leftEv(2,0) = 0.0;
  leftEv(2,1) = 0.0;
  leftEv(2,2) = 1.0;
  leftEv(2,3) = 0.0;

  leftEv(3,0) = 0.0;
  leftEv(3,1) = 0.0;
  leftEv(3,2) = 0.0;
  leftEv(3,3) = 1.0;

  eValues[0] = Vn0;
  eValues[1] = Vn1;
  eValues[2] = Vn2;
  eValues[3] = Vn3;

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @LinearAdvSys2DChar::computeEigenValuesVectors" << "\n"
  << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @LinearAdvSys2DChar::computeEigenValuesVectors" << "\n"
  << leftEv << "\n");
  CFLogDebugMax( "EigenValues @LinearAdvSys2DChar::computeEigenValuesVectors" << "\n"
  << eValues << "\n" << "\n");
}


//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setEigenVect1(RealVector& r1,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys2DPrim::setEigenVect1");
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setEigenVect2(RealVector& r2,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys2DPrim::setEigenVect2");
}
//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setEigenVect3(RealVector& r3,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys2DPrim::setEigenVect3");
} 

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setEigenVect4(RealVector& r4,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys2DPrim::setEigenVect4");
}

//////////////////////////////////////////////////////////////////////////////

CFuint LinearAdvSys2DPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  Common::SafePtr<LinearAdvSysTerm> term = getModel();
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal c0x = term->getc0x();
  const CFreal c1x = term->getc1x();
  const CFreal c2x = term->getc2x();
  const CFreal c3x = term->getc3x();

  const CFreal c0y = term->getc0y();
  const CFreal c1y = term->getc1y();
  const CFreal c2y = term->getc2y();
  const CFreal c3y = term->getc3y();

  const CFreal Vn0 = c0x*nx + c0y*ny;
  const CFreal Vn1 = c1x*nx + c1y*ny;
  const CFreal Vn2 = c2x*nx + c2y*ny;
  const CFreal Vn3 = c3x*nx + c3y*ny;

  _rightEv(0,0) = 1.0;
  _rightEv(0,1) = 0.0;
  _rightEv(0,2) = 0.0;
  _rightEv(0,3) = 0.0;

  _rightEv(1,0) = 0.0;
  _rightEv(1,1) = 1.0;
  _rightEv(1,2) = 0.0;
  _rightEv(1,3) = 0.0;

  _rightEv(2,0) = 0.0;
  _rightEv(2,1) = 0.0;
  _rightEv(2,2) = 1.0;
  _rightEv(2,3) = 0.0;

  _rightEv(3,0) = 0.0;
  _rightEv(3,1) = 0.0;
  _rightEv(3,2) = 0.0;
  _rightEv(3,3) = 1.0;



  _leftEv(0,0) = 1.0;
  _leftEv(0,1) = 0.0;
  _leftEv(0,2) = 0.0;
  _leftEv(0,3) = 0.0;

  _leftEv(1,0) = 0.0;
  _leftEv(1,1) = 1.0;
  _leftEv(1,2) = 0.0;
  _leftEv(1,3) = 0.0;

  _leftEv(2,0) = 0.0;
  _leftEv(2,1) = 0.0;
  _leftEv(2,2) = 1.0;
  _leftEv(2,3) = 0.0;

  _leftEv(3,0) = 0.0;
  _leftEv(3,1) = 0.0;
  _leftEv(3,2) = 0.0;
  _leftEv(3,3) = 1.0;

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

void LinearAdvSys2DPrim::computePhysicalData(const Framework::State& state,
  		   RealVector& data)
{

const RealVector& linearData =
    getModel()->getPhysicalData();
  
  data[LinearAdvSysTerm::u0] = state[0];
  data[LinearAdvSysTerm::u1] = state[1];
  data[LinearAdvSysTerm::u2] = state[2];
  data[LinearAdvSysTerm::u3] = state[3];
  data[LinearAdvSysTerm::C0X] = linearData[LinearAdvSysTerm::C0X];
  data[LinearAdvSysTerm::C0Y] = linearData[LinearAdvSysTerm::C0Y];
  data[LinearAdvSysTerm::C1X] = linearData[LinearAdvSysTerm::C1X];
  data[LinearAdvSysTerm::C1Y] = linearData[LinearAdvSysTerm::C1Y];
  data[LinearAdvSysTerm::C2X] = linearData[LinearAdvSysTerm::C2X];
  data[LinearAdvSysTerm::C2Y] = linearData[LinearAdvSysTerm::C2Y];
  data[LinearAdvSysTerm::C3X] = linearData[LinearAdvSysTerm::C3X];
  data[LinearAdvSysTerm::C3Y] = linearData[LinearAdvSysTerm::C3Y];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdvSys2DPrim::getSpeed(const State& state) const
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys2DPrim::getSpeed");
  return 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setDimensionalValues(const State& state, RealVector& result)
{
  result = state;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::setAdimensionalValues(const State& state, RealVector& result)
{
  result = state;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::computePerturbedStatesData
(const vector<State*>& states,
 const CFuint nbStatesInVec,
 const CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys2DPrim::computePerturbedStatesData");
}


//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys2DPrim::computeStateFromPhysicalData(const RealVector& data,
					      State& state)
{
  state[0] = data[LinearAdvSysTerm::u0];
  state[1] = data[LinearAdvSysTerm::u1];
  state[2] = data[LinearAdvSysTerm::u2];
  state[3] = data[LinearAdvSysTerm::u3];
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

