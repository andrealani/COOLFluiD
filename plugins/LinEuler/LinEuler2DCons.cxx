#include <numeric>

#include "Environment/ObjectProvider.hh"
#include "LinearizedEuler.hh"
#include "LinEuler2DCons.hh"
#include "../src/MathTools/MatrixInverter.hh"
#include "Common/NotImplementedException.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinEuler2DCons, ConvectiveVarSet, LinearizedEulerModule, 1>
lineuler2DConsProvider("LinEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

LinEuler2DCons::LinEuler2DCons(Common::SafePtr<BaseTerm> term) :
  LinEuler2DVarSet(term),
  _rightEv(4,4),
  _leftEv(4,4)
{
  ///@todo NV:Need to check if the conservative variable are that or something else

  vector<std::string> names(4);

  names[0] = "rho";
  names[1] = "rho0u";
  names[2] = "rho0v";
  names[3] = "p";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

LinEuler2DCons::~LinEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::setup()
{
  LinEuler2DVarSet::setup();
  setConstJacob();
}
//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::setConstJacob()
{

 vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[XX](0,1) = 1.0;
  (*jacobians)[XX](0,2) = 0.0;
  (*jacobians)[XX](0,3) = 0.0;

  (*jacobians)[XX](1,0) = 0.0;
  (*jacobians)[XX](1,2) = 0.0;
  (*jacobians)[XX](1,3) = 1.0;

  (*jacobians)[XX](2,0) = 0.0;
  (*jacobians)[XX](2,1) = 0.0;
  (*jacobians)[XX](2,3) = 0.0;

  (*jacobians)[XX](3,0) = 0.0;
  (*jacobians)[XX](3,2) = 0.0;

  (*jacobians)[YY](0,1) = 0.0;
  (*jacobians)[YY](0,2) = 1.0;
  (*jacobians)[YY](0,3) = 0.0;

  (*jacobians)[YY](1,0) = 0.0;
  (*jacobians)[YY](1,2) = 0.0;
  (*jacobians)[YY](1,3) = 0.0;

  (*jacobians)[YY](2,0) = 0.0;
  (*jacobians)[YY](2,1) = 0.0;
  (*jacobians)[YY](2,3) = 1.0;

  (*jacobians)[YY](3,0) = 0.0;
  (*jacobians)[YY](3,1) = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::computeProjectedJacobian(const RealVector& normal,
            RealMatrix& jacob)
{
  RealVector& linearData = _model->getPhysicalData();


  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal c     = linearData[LinEulerTerm::c];
  
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal Vn = U0*nx+V0*ny;

  jacob(0,1) = 1.0*nx;
  jacob(0,2) = 1.0*ny;
  jacob(0,3) = 0.0;
  jacob(1,0) = 0.0;
  jacob(1,2) = 0.0;
  jacob(1,3) = 1.0*nx;
  jacob(2,0) = 0.0;
  jacob(2,1) = 0.0;
  jacob(2,3) = 1.0*ny;
  jacob(3,0) = 0.0;

  jacob(0,0) = Vn;
  jacob(1,1) = Vn;
  jacob(2,2) = Vn;
  jacob(3,3) = Vn;
  jacob(3,1) = c*c*nx;
  jacob(3,2) = c*c*ny;

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::computeJacobians()
{
  vector<RealMatrix>* const jacobians =     PhysicalModelStack::getActive()->getImplementor()->getJacobians();

//  Common::SafePtr<LinEulerTerm> linearData = getModel();
  RealVector& linearData = _model->getPhysicalData();

  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal c     = linearData[LinEulerTerm::c];

/*
   CF_DEBUG_OBJ(U0);
   CF_DEBUG_OBJ(V0);
   CF_DEBUG_OBJ(gamma);
   CF_DEBUG_OBJ(c);
*/


  (*jacobians)[XX](0,0) = U0;
  (*jacobians)[XX](1,1) = U0;
  (*jacobians)[XX](2,2) = U0;
  (*jacobians)[XX](3,1) = c*c;
  (*jacobians)[XX](3,3) = U0;
  (*jacobians)[YY](0,0) = V0;
  (*jacobians)[YY](1,1) = V0;
  (*jacobians)[YY](2,2) = V0;
  (*jacobians)[YY](3,2) = c*c;
  (*jacobians)[YY](3,3) = V0;

}


//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::computeEigenValuesVectors(RealMatrix& rightEv,
          RealMatrix& leftEv,
          RealVector& eValues,
          const RealVector& normal) {

/* the implementation waiting for normalized normal vectors!!! */

  RealVector& linearData = _model->getPhysicalData();

  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal Vn = U0*nx+V0*ny;

  const CFreal inv_c = 1./c;
  const CFreal inv_c2  = inv_c/c;
/*
   CF_DEBUG_OBJ(U0);
   CF_DEBUG_OBJ(V0);
   CF_DEBUG_OBJ(c);
*/

  rightEv(0,0) = 1.0;
  rightEv(0,1) = 0.0;
  rightEv(0,2) = 0.5*inv_c;
  rightEv(0,3) = 0.5*inv_c;

  rightEv(1,0) = 0.0;
  rightEv(1,1) = ny;
  rightEv(1,2) = 0.5*nx;
  rightEv(1,3) = -0.5*nx;

  rightEv(2,0) = 0.0;
  rightEv(2,1) = -nx;
  rightEv(2,2) = 0.5*ny;
  rightEv(2,3) = -0.5*ny;

  rightEv(3,0) = 0.0;
  rightEv(3,1) = 0.0;
  rightEv(3,2) = 0.5*c;
  rightEv(3,3) = 0.5*c;

 leftEv(0,0) = 1.0;
 leftEv(0,1) = 0.0;
 leftEv(0,2) = 0.0;
 leftEv(0,3) = -inv_c2;

 leftEv(1,0) = 0.0;
 leftEv(1,1) = ny;
 leftEv(1,2) = -nx;
 leftEv(1,3) = 0.0;

 leftEv(2,0) = 0.0;
 leftEv(2,1) = nx;
 leftEv(2,2) = ny;
 leftEv(2,3) = inv_c;

 leftEv(3,0) = 0.0;
 leftEv(3,1) = -nx;
 leftEv(3,2) = -ny;
 leftEv(3,3) = inv_c;

 eValues[0] = Vn;
 eValues[1] = Vn;
 eValues[2] = Vn+c;
 eValues[3] = Vn-c;

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @LinEuler2DChar::setEigenValuesVectors" << "\n"
  << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @LinEuler2DChar::setEigenValuesVectors" << "\n"
  << leftEv << "\n");
  CFLogDebugMax( "EigenValues @LinEuler2DChar::setEigenValuesVectors" << "\n"
  << eValues << "\n" << "\n");
}


//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::setEigenVect1(RealVector& r1,
        State& state,
        const RealVector& normal)
{

  r1[0] = 1.0;
  r1[1] = 0.0;
  r1[2] = 0.0;
  r1[3] = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::setEigenVect2(RealVector& r2,
        State& state,
        const RealVector& normal)
{

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  r2[0] = 0.0;
  r2[1] = ny;
  r2[2] = -nx;
  r2[3] = 0.0;

}
//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::setEigenVect3(RealVector& r3,
        State& state,
        const RealVector& normal)
{

  RealVector& linearData = _model->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  r3[0] = 0.5/c;
  r3[1] = 0.5*nx;
  r3[2] = 0.5*ny;
  r3[3] = 0.5*c;

}

//////////////////////////////////////////////////////////////////////////////
void LinEuler2DCons::setEigenVect4(RealVector& r4,
        State& state,
        const RealVector& normal)
{

  RealVector& linearData = _model->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  r4[0] = 0.5/c;
  r4[1] = -0.5*nx;
  r4[2] = -0.5*ny;
  r4[3] = 0.5*c;

}

//////////////////////////////////////////////////////////////////////////////

CFuint LinEuler2DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal Vn = U0*nx+V0*ny;

/*
   CF_DEBUG_OBJ(U0);
   CF_DEBUG_OBJ(V0);
   CF_DEBUG_OBJ(c);
*/

  const CFreal inv_c = 1./c;
  const CFreal inv_c2  = inv_c/c;

  _rightEv(0,0) = 1.0;
  _rightEv(0,1) = 0.0;
  _rightEv(0,2) = 0.5*inv_c;
  _rightEv(0,3) = 0.5*inv_c;

  _rightEv(1,0) = 0.0;
  _rightEv(1,1) = ny;
  _rightEv(1,2) = 0.5*nx;
  _rightEv(1,3) = -0.5*nx;

  _rightEv(2,0) = 0.0;
  _rightEv(2,1) = -nx;
  _rightEv(2,2) = 0.5*ny;
  _rightEv(2,3) = -0.5*ny;

  _rightEv(3,0) = 0.0;
  _rightEv(3,1) = 0.0;
  _rightEv(3,2) = 0.5*c;
  _rightEv(3,3) = 0.5*c;

  _leftEv(0,0) = 1.0;
  _leftEv(0,1) = 0.0;
  _leftEv(0,2) = 0.0;
  _leftEv(0,3) = -inv_c2;

  _leftEv(1,0) = 0.0;
  _leftEv(1,1) = ny;
  _leftEv(1,2) = -nx;
  _leftEv(1,3) = 0.0;

  _leftEv(2,0) = 0.0;
  _leftEv(2,1) = nx;
  _leftEv(2,2) = ny;
  _leftEv(2,3) = inv_c;

  _leftEv(3,0) = 0.0;
  _leftEv(3,1) = -nx;
  _leftEv(3,2) = -ny;
  _leftEv(3,3) = inv_c;

  eValues[0] = Vn;
  eValues[1] = Vn;
  eValues[2] = Vn+c;
  eValues[3] = Vn-c;


  for (CFuint iEq = 0; iEq < eValues.size(); ++iEq) {
      _eValuesP[iEq] = max(0.,eValues[iEq]);
      _eValuesM[iEq] = min(0.,eValues[iEq]);
  }

  jacobPlus = _rightEv*(_eValuesP*_leftEv);
  jacobMin  = _rightEv*(_eValuesM*_leftEv);

  CFLogDebugMax( "RightEigenvectors @LinEuler2DCons::splitJacob" << "\n"
  << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @LinEuler2DCons::splitJacob" << "\n"
  << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @LinEuler2DCons::splitJacob" << "\n"
  << eValues << "\n" << "\n");

}//End of LinEuler2DCons::splitJacob

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::computePhysicalData(const State& state, RealVector& data)
{
  RealVector& linearData = _model->getPhysicalData();
  
  const CFreal rho0  = linearData[LinEulerTerm::rho0];
  const CFreal rho  = state[0];
  const CFreal u = state[1]/rho0;
  const CFreal v = state[2]/rho0;
  const CFreal p = state[3];

  data[LinEulerTerm::rho] = rho;
  data[LinEulerTerm::u]   = u;
  data[LinEulerTerm::v]   = v;
  data[LinEulerTerm::p]   = p;

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::computeStateFromPhysicalData(const RealVector& data, State& state)
{
  RealVector& linearData = _model->getPhysicalData();
  const CFreal rho0  = linearData[LinEulerTerm::rho0];

  state[0] = data[LinEulerTerm::rho];
  state[1] = rho0*data[LinEulerTerm::u];
  state[2] = rho0*data[LinEulerTerm::v];
  state[3] = data[LinEulerTerm::p];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinEuler2DCons::getSpeed(const State& state) const
{

  RealVector& linearData = _model->getPhysicalData();
  const CFreal rho0  = linearData[LinEulerTerm::rho0];

  const CFreal u = state[1]/rho0;
  const CFreal v = state[2]/rho0;
  return sqrt(u*u + v*v);

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  result[0] = state[0]*1.0;
  result[1] = state[1]*1.0;
  result[2] = state[2]*1.0;
  result[3] = state[3]*1.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  result[0] = state[0]*1.0;
  result[1] = state[1]*1.0;
  result[2] = state[2]*1.0;
  result[3] = state[3]*1.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DCons::computePerturbedStatesData
(const vector<State*>& states,
 const CFuint nbStatesInVec,
 const CFuint iVar)
{
throw Common::NotImplementedException(FromHere(),"LinEuler2DCons::computePerturbedStatesData");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
