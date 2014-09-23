#include <numeric>

#include "Environment/ObjectProvider.hh"
#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/LinEuler3DCons.hh"
#include "../src/MathTools/MatrixInverter.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinEuler3DCons, ConvectiveVarSet, LinearizedEulerModule, 1>
lineuler3DConsProvider("LinEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

LinEuler3DCons::LinEuler3DCons(SafePtr<BaseTerm> term) :
  LinEuler3DVarSet(term),
  _rightEv(5,5),
  _leftEv(5,5)
{
  vector<std::string> names(5);

  names[0] = "rho";
  names[1] = "rho0u";
  names[2] = "rho0v";
  names[3] = "rho0w";
  names[4] = "p";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

LinEuler3DCons::~LinEuler3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setup()
{
  LinEuler3DVarSet::setup();
  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setConstJacob()
{

 vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[XX](0,1) = 1.0;
  (*jacobians)[XX](0,2) = 0.0;
  (*jacobians)[XX](0,3) = 0.0;
  (*jacobians)[XX](0,4) = 0.0;

  (*jacobians)[XX](1,0) = 0.0;
  (*jacobians)[XX](1,2) = 0.0;
  (*jacobians)[XX](1,3) = 0.0;
  (*jacobians)[XX](1,4) = 1.0;

  (*jacobians)[XX](2,0) = 0.0;
  (*jacobians)[XX](2,1) = 0.0;
  (*jacobians)[XX](2,3) = 0.0;
  (*jacobians)[XX](2,4) = 0.0;

  (*jacobians)[XX](3,0) = 0.0;
  (*jacobians)[XX](3,1) = 0.0;
  (*jacobians)[XX](3,2) = 0.0;
  (*jacobians)[XX](3,4) = 0.0;

  (*jacobians)[XX](4,0) = 0.0;
  (*jacobians)[XX](4,2) = 0.0;
  (*jacobians)[XX](4,3) = 0.0;


  (*jacobians)[YY](0,1) = 0.0;
  (*jacobians)[YY](0,2) = 1.0;
  (*jacobians)[YY](0,3) = 0.0;
  (*jacobians)[YY](0,4) = 0.0;

  (*jacobians)[YY](1,0) = 0.0;
  (*jacobians)[YY](1,2) = 0.0;
  (*jacobians)[YY](1,3) = 0.0;
  (*jacobians)[YY](1,4) = 0.0;

  (*jacobians)[YY](2,0) = 0.0;
  (*jacobians)[YY](2,1) = 0.0;
  (*jacobians)[YY](2,3) = 0.0;
  (*jacobians)[YY](2,4) = 1.0;

  (*jacobians)[YY](3,0) = 0.0;
  (*jacobians)[YY](3,1) = 0.0;
  (*jacobians)[YY](3,2) = 0.0;
  (*jacobians)[YY](3,4) = 0.0;

  (*jacobians)[YY](4,0) = 0.0;
  (*jacobians)[YY](4,1) = 0.0;
  (*jacobians)[YY](4,3) = 0.0;


  (*jacobians)[ZZ](0,1) = 0.0;
  (*jacobians)[ZZ](0,2) = 0.0;
  (*jacobians)[ZZ](0,3) = 1.0;
  (*jacobians)[ZZ](0,4) = 0.0;

  (*jacobians)[ZZ](1,0) = 0.0;
  (*jacobians)[ZZ](1,2) = 0.0;
  (*jacobians)[ZZ](1,3) = 0.0;
  (*jacobians)[ZZ](1,4) = 0.0;

  (*jacobians)[ZZ](2,0) = 0.0;
  (*jacobians)[ZZ](2,1) = 0.0;
  (*jacobians)[ZZ](2,3) = 0.0;
  (*jacobians)[ZZ](2,4) = 0.0;

  (*jacobians)[ZZ](3,0) = 0.0;
  (*jacobians)[ZZ](3,1) = 0.0;
  (*jacobians)[ZZ](3,2) = 0.0;
  (*jacobians)[ZZ](3,4) = 1.0;

  (*jacobians)[ZZ](4,0) = 0.0;
  (*jacobians)[ZZ](4,1) = 0.0;
  (*jacobians)[ZZ](4,2) = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::computeProjectedJacobian(const RealVector& normal,
            RealMatrix& jacob)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal W0    = linearData[LinEulerTerm::W0];
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal Vn = U0*nx+V0*ny+W0*nz;


  jacob(0,0) = Vn;
  jacob(0,1) = nx;
  jacob(0,2) = ny;
  jacob(0,3) = nz;
  jacob(0,4) = 0.0;

  jacob(1,0) = 0.0;
  jacob(1,1) = Vn;
  jacob(1,2) = 0.0;
  jacob(1,3) = 0.0;
  jacob(1,4) = nx;

  jacob(2,0) = 0.0;
  jacob(2,1) = 0.0;
  jacob(2,2) = Vn;
  jacob(2,3) = 0.0;
  jacob(2,4) = ny;

  jacob(3,0) = 0.0;
  jacob(3,1) = 0.0;
  jacob(3,2) = 0.0;
  jacob(3,3) = Vn;
  jacob(3,4) = nz;

  jacob(4,0) = 0.0;
  jacob(4,1) = c*c*nx;
  jacob(4,2) = c*c*ny;
  jacob(4,3) = c*c*nz;
  jacob(4,4) = Vn;

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::computeJacobians()
{
  vector<RealMatrix>* const jacobians =     PhysicalModelStack::getActive()->getImplementor()->getJacobians();

RealVector& linearData = _model->getPhysicalData();

  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal W0    = linearData[LinEulerTerm::W0];
  const CFreal c     = linearData[LinEulerTerm::c];

  (*jacobians)[XX](0,0) = U0;
  (*jacobians)[XX](1,1) = U0;
  (*jacobians)[XX](2,2) = U0;
  (*jacobians)[XX](3,3) = U0;
  (*jacobians)[XX](4,4) = U0;

  (*jacobians)[XX](4,1) = c*c;

  (*jacobians)[YY](0,0) = V0;
  (*jacobians)[YY](1,1) = V0;
  (*jacobians)[YY](2,2) = V0;
  (*jacobians)[YY](3,3) = V0;
  (*jacobians)[YY](4,4) = V0;

  (*jacobians)[YY](4,2) = c*c;

  (*jacobians)[ZZ](0,0) = W0;
  (*jacobians)[ZZ](1,1) = W0;
  (*jacobians)[ZZ](2,2) = W0;
  (*jacobians)[ZZ](3,3) = W0;
  (*jacobians)[ZZ](4,4) = W0;

  (*jacobians)[ZZ](4,3) = c*c;


}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                               RealMatrix& leftEv,
                                               RealVector& eValues,
                                               const RealVector& normal)
{
/* the implementation waiting for normalized normal vectors!!! */

  RealVector& linearData = _model->getPhysicalData();

  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal W0    = linearData[LinEulerTerm::W0];
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal Vn = U0*nx+V0*ny+W0*nz;

  const CFreal inv_c = 1./c;
  const CFreal inv_c2  = inv_c/c;

  rightEv(0,0) = nx;
  rightEv(0,1) = ny;
  rightEv(0,2) = nz;
  rightEv(0,3) = 0.5*inv_c;
  rightEv(0,4) = 0.5*inv_c;

  rightEv(1,0) = 0.0;
  rightEv(1,1) = -nz;
  rightEv(1,2) = ny;
  rightEv(1,3) = 0.5*nx;
  rightEv(1,4) = -0.5*nx;

  rightEv(2,0) = nz;
  rightEv(2,1) = 0.0;
  rightEv(2,2) = -nx;
  rightEv(2,3) = 0.5*ny;
  rightEv(2,4) = -0.5*ny;

  rightEv(3,0) = -ny;
  rightEv(3,1) = nx;
  rightEv(3,2) = 0.0;
  rightEv(3,3) = 0.5*nz;
  rightEv(3,4) = -0.5*nz;

  rightEv(4,0) = 0.0;
  rightEv(4,1) = 0.0;
  rightEv(4,2) = 0.0;
  rightEv(4,3) = 0.5*c;
  rightEv(4,4) = 0.5*c;


 leftEv(0,0) = nx;
 leftEv(0,1) = 0.0;
 leftEv(0,2) = nz;
 leftEv(0,3) = -ny;
 leftEv(0,4) = -nx*inv_c2;

 leftEv(1,0) = ny;
 leftEv(1,1) = -nz;
 leftEv(1,2) = 0.0;
 leftEv(1,3) = nx;
 leftEv(1,4) = -ny*inv_c2;

 leftEv(2,0) = nz;
 leftEv(2,1) = ny;
 leftEv(2,2) = -nx;
 leftEv(2,3) = 0.0;
 leftEv(2,4) = -nz*inv_c2;

 leftEv(3,0) = 0.0;
 leftEv(3,1) = nx;
 leftEv(3,2) = ny;
 leftEv(3,3) = nz;
 leftEv(3,4) = inv_c;

 leftEv(4,0) = 0.0;
 leftEv(4,1) = -nx;
 leftEv(4,2) = -ny;
 leftEv(4,3) = -nz;
 leftEv(4,4) = inv_c;

 eValues[0] = Vn;
 eValues[1] = Vn;
 eValues[2] = Vn;
 eValues[3] = Vn+c;
 eValues[4] = Vn-c;

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @LinEuler2DChar::setEigenValuesVectors" << "\n"
  << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @LinEuler2DChar::setEigenValuesVectors" << "\n"
  << leftEv << "\n");
  CFLogDebugMax( "EigenValues @LinEuler2DChar::setEigenValuesVectors" << "\n"
  << eValues << "\n" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setEigenVect1(RealVector& r1,
                                   State& state,
                                   const RealVector& normal)
{

  const CFreal nz = normal[ZZ];
  const CFreal ny = normal[YY];
  const CFreal nx = normal[XX];

  r1[0] = nx;
  r1[1] = 0.0;
  r1[2] = nz;
  r1[3] = -ny;
  r1[4] = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setEigenVect2(RealVector& r2,
                                   State& state,
                                   const RealVector& normal)
{
  const CFreal nz = normal[ZZ];
  const CFreal ny = normal[YY];
  const CFreal nx = normal[XX];

  r2[0] = ny;
  r2[1] = -nz;
  r2[2] = 0.0;
  r2[3] = nx;
  r2[4] = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setEigenVect3(RealVector& r3,
                                   State& state,
                                   const RealVector& normal)
{
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];

  r3[0] = nz;
  r3[1] = ny;
  r3[2] = -nx;
  r3[3] = 0.0;
  r3[4] = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setEigenVect4(RealVector& r4,
                                   State& state,
                                   const RealVector& normal)
{
  RealVector& linearData = _model->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];

  r4[0] = 0.5/c;
  r4[1] = 0.5*nx;
  r4[2] = 0.5*ny;
  r4[3] = 0.5*nz;
  r4[4] = 0.5*c;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setEigenVect5(RealVector& r5,
                                   State& state,
                                   const RealVector& normal)
{
  RealVector& linearData = _model->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];

  r5[0] = 0.5/c;
  r5[1] = -0.5*nx;
  r5[2] = -0.5*ny;
  r5[3] = -0.5*nz;
  r5[4] = -0.5*c;
}

//////////////////////////////////////////////////////////////////////////////

CFuint LinEuler3DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::splitJacobian(RealMatrix& jacobPlus,
                                   RealMatrix& jacobMin,
                                   RealVector& eValues,
                                   const RealVector& normal)
{
 RealVector& linearData = _model->getPhysicalData();

  const CFreal U0    = linearData[LinEulerTerm::U0];
  const CFreal V0    = linearData[LinEulerTerm::V0];
  const CFreal W0    = linearData[LinEulerTerm::W0];
  const CFreal c     = linearData[LinEulerTerm::c];


  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal Vn = U0*nx+V0*ny+W0*nz;

  const CFreal inv_c = 1./c;
  const CFreal inv_c2  = inv_c/c;

_rightEv(0,0) = nx;
_rightEv(0,1) = ny;
_rightEv(0,2) = nz;
_rightEv(0,3) = 0.5*inv_c;
_rightEv(0,4) = 0.5*inv_c;

_rightEv(1,0) = 0.0;
_rightEv(1,1) = -nz;
_rightEv(1,2) = ny;
_rightEv(1,3) = 0.5*nx;
_rightEv(1,4) = -0.5*nx;

_rightEv(2,0) = nz;
_rightEv(2,1) = 0.0;
_rightEv(2,2) = -nx;
_rightEv(2,3) = 0.5*ny;
_rightEv(2,4) = -0.5*ny;

_rightEv(3,0) = -ny;
_rightEv(3,1) = nx;
_rightEv(3,2) = 0.0;
_rightEv(3,3) = 0.5*nz;
_rightEv(3,4) = -0.5*nz;

_rightEv(4,0) = 0.0;
_rightEv(4,1) = 0.0;
_rightEv(4,2) = 0.0;
_rightEv(4,3) = 0.5*c;
_rightEv(4,4) = 0.5*c;


_leftEv(0,0) = nx;
_leftEv(0,1) = 0.0;
_leftEv(0,2) = nz;
_leftEv(0,3) = -ny;
_leftEv(0,4) = -nx*inv_c2;

_leftEv(1,0) = ny;
_leftEv(1,1) = -nz;
_leftEv(1,2) = 0.0;
_leftEv(1,3) = nx;
_leftEv(1,4) = -ny*inv_c2;

_leftEv(2,0) = nz;
_leftEv(2,1) = ny;
_leftEv(2,2) = -nx;
_leftEv(2,3) = 0.0;
_leftEv(2,4) = -nz*inv_c2;

_leftEv(3,0) = 0.0;
_leftEv(3,1) = nx;
_leftEv(3,2) = ny;
_leftEv(3,3) = nz;
_leftEv(3,4) = inv_c;

_leftEv(4,0) = 0.0;
_leftEv(4,1) = -nx;
_leftEv(4,2) = -ny;
_leftEv(4,3) = -nz;
_leftEv(4,4) = inv_c;

 eValues[0] = Vn;
 eValues[1] = Vn;
 eValues[2] = Vn;
 eValues[3] = Vn+c;
 eValues[4] = Vn-c;


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

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::computePhysicalData(const State& state,
                                         RealVector& data)
{
  // get local (reconstructed) mean flow

  RealVector& linearData = _model->getPhysicalData();

  const CFreal rho0 = linearData[LinEulerTerm::rho0];
  const CFreal U0   = linearData[LinEulerTerm::U0];
  const CFreal V0   = linearData[LinEulerTerm::V0];
  const CFreal W0   = linearData[LinEulerTerm::W0];
  const CFreal P0   = linearData[LinEulerTerm::P0];
  const CFreal gamma = _model->getgamma();

  // set mean flow in data
  data[LinEulerTerm::GAMMA] = gamma;
  data[LinEulerTerm::rho0]  = rho0;
  data[LinEulerTerm::U0]    = U0;
  data[LinEulerTerm::V0]    = V0;
  data[LinEulerTerm::W0]    = W0;
  data[LinEulerTerm::P0]    = P0;
  data[LinEulerTerm::c]     = sqrt(gamma*P0/rho0);

  // set perturbations
  data[LinEulerTerm::rho] = state[0];
  data[LinEulerTerm::u]   = state[1]/rho0;
  data[LinEulerTerm::v]   = state[2]/rho0;
  data[LinEulerTerm::w]   = state[3]/rho0;
  data[LinEulerTerm::p]   = state[4];
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::computeStateFromPhysicalData(const RealVector& data,
                                                  State& state)
{

  RealVector& linearData = _model->getPhysicalData();


  cf_assert(_localMeanFlow != CFNULL);
  const CFreal rho0 = linearData[LinEulerTerm::rho0];

  state[0] = data[LinEulerTerm::rho];
  state[1] = rho0*data[LinEulerTerm::u];
  state[2] = rho0*data[LinEulerTerm::v];
  state[3] = rho0*data[LinEulerTerm::w];
  state[4] = data[LinEulerTerm::p];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinEuler3DCons::getSpeed(const State& state) const
{
  RealVector& linearData = _model->getPhysicalData();
  const CFreal rho0  = linearData[LinEulerTerm::rho0];

  const CFreal u = state[1]/rho0;
  const CFreal v = state[2]/rho0;
  const CFreal w = state[3]/rho0;
  return sqrt(u*u + v*v + w*w);
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  result[0] = state[0]*1.0;
  result[1] = state[1]*1.0;
  result[2] = state[2]*1.0;
  result[3] = state[3]*1.0;
  result[4] = state[4]*1.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::setAdimensionalValues(const State& state,
                                           RealVector& result)
{
  result[0] = state[0]*1.0;
  result[1] = state[1]*1.0;
  result[2] = state[2]*1.0;
  result[3] = state[3]*1.0;
  result[4] = state[4]*1.0;
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DCons::computePerturbedStatesData(const vector<State*>& states,
                                                const CFuint nbStatesInVec,
                                                const CFuint iVar)
{
  throw NotImplementedException(FromHere(),"LinEuler3DCons::computePerturbedStatesData");
}

// //////////////////////////////////////////////////////////////////////////////
// 
// void LinEuler3DCons::computeFlux(const RealVector& pdata, const RealVector& normals)
// {
//  throw Common::NotImplementedException(FromHere(), "LinEuler3DCons::computeFlux()");
// }

// //////////////////////////////////////////////////////////////////////////////
// 
// void LinEuler3DCons::computeStateFlux(const RealVector& vars)
// {
//  throw Common::NotImplementedException(FromHere(), "LinEuler3DCons::computeStateFlux()");
// }

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

