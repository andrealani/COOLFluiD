// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "LinearAdvSys/LinearAdvSys.hh"
#include "LinearAdvSys/LinearAdvSys3DPrim.hh"
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

Environment::ObjectProvider<LinearAdvSys3DPrim, ConvectiveVarSet, LinearAdvSysModule, 1>
linearAdvSys3DPrimProvider("LinearAdvSys3DPrim");

//////////////////////////////////////////////////////////////////////////////

LinearAdvSys3DPrim::LinearAdvSys3DPrim(Common::SafePtr<BaseTerm> term) :
  LinearAdvSys3DVarSet(term),
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

LinearAdvSys3DPrim::~LinearAdvSys3DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setup()
{
  LinearAdvSys3DVarSet::setup();
  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setConstJacob()
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

  const CFreal c0z = term->getc0z();
  const CFreal c1z = term->getc1z();
  const CFreal c2z = term->getc2z();
  const CFreal c3z = term->getc3z();

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

  (*jacobians)[ZZ](0,0) = c0z;
  (*jacobians)[ZZ](0,1) = 0.0;
  (*jacobians)[ZZ](0,2) = 0.0;
  (*jacobians)[ZZ](0,3) = 0.0;

  (*jacobians)[ZZ](1,0) = 0.0;
  (*jacobians)[ZZ](1,1) = c1z;
  (*jacobians)[ZZ](1,2) = 0.0;
  (*jacobians)[ZZ](1,3) = 0.0;

  (*jacobians)[ZZ](2,0) = 0.0;
  (*jacobians)[ZZ](2,1) = 0.0;
  (*jacobians)[ZZ](2,2) = c2z;
  (*jacobians)[ZZ](2,3) = 0.0;

  (*jacobians)[ZZ](3,0) = 0.0;
  (*jacobians)[ZZ](3,1) = 0.0;
  (*jacobians)[ZZ](3,2) = 0.0;
  (*jacobians)[ZZ](3,3) = c3z;

}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob)
{
  /// @note TK: Because of the Jacobian is given by main flow it is constnt
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::computeJacobians()
{
  /// @note TK: Because of the Jacobian is given by main flow it is constnt
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
          RealMatrix& leftEv,
          RealVector& eValues,
          const RealVector& normal)
{
  Common::SafePtr<LinearAdvSysTerm> term = getModel();
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];

  const CFreal c0x = term->getc0x();
  const CFreal c1x = term->getc1x();
  const CFreal c2x = term->getc2x();
  const CFreal c3x = term->getc3x();

  const CFreal c0y = term->getc0y();
  const CFreal c1y = term->getc1y();
  const CFreal c2y = term->getc2y();
  const CFreal c3y = term->getc3y();

  const CFreal c0z = term->getc0z();
  const CFreal c1z = term->getc1z();
  const CFreal c2z = term->getc2z();
  const CFreal c3z = term->getc3z();

  const CFreal Vn0 = c0x*nx + c0y*ny + c0z*nz;
  const CFreal Vn1 = c1x*nx + c1y*ny + c1z*nz;
  const CFreal Vn2 = c2x*nx + c2y*ny + c2z*nz;
  const CFreal Vn3 = c3x*nx + c3y*ny + c3z*nz;

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
  CFLogDebugMax( "RightEigenvectors @LinearAdvSys3DChar::computeEigenValuesVectors" << "\n"
  << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @LinearAdvSys3DChar::computeEigenValuesVectors" << "\n"
  << leftEv << "\n");
  CFLogDebugMax( "EigenValues @LinearAdvSys3DChar::computeEigenValuesVectors" << "\n"
  << eValues << "\n" << "\n");
}


//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setEigenVect1(RealVector& r1,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys3DPrim::setEigenVect1");
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setEigenVect2(RealVector& r2,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys3DPrim::setEigenVect2");
}
//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setEigenVect3(RealVector& r3,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys3DPrim::setEigenVect3");
} 

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setEigenVect4(RealVector& r4,
        State& state,
        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys3DPrim::setEigenVect4");
}

//////////////////////////////////////////////////////////////////////////////

CFuint LinearAdvSys3DPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  Common::SafePtr<LinearAdvSysTerm> term = getModel();
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];

  const CFreal c0x = term->getc0x();
  const CFreal c1x = term->getc1x();
  const CFreal c2x = term->getc2x();
  const CFreal c3x = term->getc3x();

  const CFreal c0y = term->getc0y();
  const CFreal c1y = term->getc1y();
  const CFreal c2y = term->getc2y();
  const CFreal c3y = term->getc3y();

  const CFreal c0z = term->getc0z();
  const CFreal c1z = term->getc1z();
  const CFreal c2z = term->getc2z();
  const CFreal c3z = term->getc3z();

  const CFreal Vn0 = c0x*nx + c0y*ny + c0z*nz;
  const CFreal Vn1 = c1x*nx + c1y*ny + c1z*nz;
  const CFreal Vn2 = c2x*nx + c2y*ny + c2z*nz;
  const CFreal Vn3 = c3x*nx + c3y*ny + c3z*nz;
                                            
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

void LinearAdvSys3DPrim::computePhysicalData(const Framework::State& state,
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
  data[LinearAdvSysTerm::C0Z] = linearData[LinearAdvSysTerm::C0Z];
  data[LinearAdvSysTerm::C1X] = linearData[LinearAdvSysTerm::C1X];
  data[LinearAdvSysTerm::C1Y] = linearData[LinearAdvSysTerm::C1Y];
  data[LinearAdvSysTerm::C1Z] = linearData[LinearAdvSysTerm::C1Z];
  data[LinearAdvSysTerm::C2X] = linearData[LinearAdvSysTerm::C2X];
  data[LinearAdvSysTerm::C2Y] = linearData[LinearAdvSysTerm::C2Y];
  data[LinearAdvSysTerm::C2Z] = linearData[LinearAdvSysTerm::C2Z];
  data[LinearAdvSysTerm::C3X] = linearData[LinearAdvSysTerm::C3X];
  data[LinearAdvSysTerm::C3Y] = linearData[LinearAdvSysTerm::C3Y];
  data[LinearAdvSysTerm::C3Z] = linearData[LinearAdvSysTerm::C3Z];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinearAdvSys3DPrim::getSpeed(const State& state) const
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys3DPrim::getSpeed");
  return 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setDimensionalValues(const State& state, RealVector& result)
{
  result = state;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::setAdimensionalValues(const State& state, RealVector& result)
{
  result = state;
}

//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::computePerturbedStatesData
(const vector<State*>& states,
 const CFuint nbStatesInVec,
 const CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"LinearAdvSys3DPrim::computePerturbedStatesData");
}


//////////////////////////////////////////////////////////////////////////////

void LinearAdvSys3DPrim::computeStateFromPhysicalData(const RealVector& data,
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

