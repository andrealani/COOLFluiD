#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler2DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

LinEuler2DVarSet::LinEuler2DVarSet(Common::SafePtr<BaseTerm> term) :
  LinEulerVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

LinEuler2DVarSet::~LinEuler2DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DVarSet::setup()
{
  LinEulerVarSet::setup();

  // set EquationSetData
  LinEuler2DVarSet::getEqSetData().resize(1);
  LinEuler2DVarSet::getEqSetData()[0].setup(0,0,4);
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinEuler2DVarSet::getMaxEigenValue(const RealVector& pdata,
               const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return un + pdata[LinEulerTerm::c];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinEuler2DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
            const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return fabs(un) + pdata[LinEulerTerm::c];
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DVarSet::computeEigenValues(const RealVector& pdata,
               const RealVector& normal, RealVector& result)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();

  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID())
  {
    const vector<CFuint>& varIDs = LinEuler2DVarSet::getEqSetData()[0].getEqSetVarIDs();
    const CFreal un = getNormalSpeed(pdata, normal);
    if (varIDs.size() == 4) {

      result[varIDs[0]] = un;
      result[varIDs[1]] = un;
      result[varIDs[2]] = un + pdata[LinEulerTerm::c];
      result[varIDs[3]] = un - pdata[LinEulerTerm::c];
    }
    else if (varIDs.size() == 3) {

      result[varIDs[0]] = un;
      result[varIDs[1]] = un + pdata[LinEulerTerm::c];
      result[varIDs[2]] = un - pdata[LinEulerTerm::c];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DVarSet::computeFlux (const RealVector& pdata, const RealVector& normals)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();

  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID())
  {

//     const CFreal gamma = _model->getgamma();
//     const CFreal rho0  = pdata[LinEulerTerm::c];
//     const CFreal P0    = pdata[LinEulerTerm::P0];
//     const CFreal U0    = pdata[LinEulerTerm::U0];
//     const CFreal V0    = pdata[LinEulerTerm::V0];
    
    
    RealVector& linearData = _model->getPhysicalData();
    const CFreal gamma = _model->getgamma();
    const CFreal rho0  = linearData[LinEulerTerm::rho0];;
    const CFreal P0    = linearData[LinEulerTerm::P0];;
    const CFreal U0    = linearData[LinEulerTerm::U0];
    const CFreal V0    = linearData[LinEulerTerm::V0];
    
    

    const CFreal nx = normals[XX];
    const CFreal ny = normals[YY];
    const CFreal V0n = U0*nx+V0*ny;
    const CFreal u = pdata[LinEulerTerm::u];
    const CFreal v = pdata[LinEulerTerm::v];
    const CFreal un = u*nx + v*ny;
    const CFreal rho = pdata[LinEulerTerm::rho];
    const CFreal p = pdata[LinEulerTerm::p];

 EquationSetData& eqSet = LinEuler2DVarSet::getEqSetData()[0];
 const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();


/*
   CF_DEBUG_OBJ(U0);
   CF_DEBUG_OBJ(V0);
   CF_DEBUG_OBJ(gamma);
   CF_DEBUG_OBJ(rho0);
   CF_DEBUG_OBJ(P0);
*/
 

 if (varIDs.size() == 4) {
     _fluxArray[varIDs[0]] = V0n*rho+un*rho0;
     _fluxArray[varIDs[1]] = rho0*V0n*u+p*nx;
     _fluxArray[varIDs[2]] = rho0*V0n*v+p*ny;
     _fluxArray[varIDs[3]] = V0n*p+un*gamma*P0;
    }
    else if (varIDs.size() == 3) {
      CF_DEBUG_EXIT;
     _fluxArray[varIDs[0]] = rho0*V0n*u+p*nx;
     _fluxArray[varIDs[1]] = rho0*V0n*v+p*ny;
     _fluxArray[varIDs[2]] = V0n*p+un*gamma*P0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler2DVarSet::computeStateFlux (const RealVector& pdata)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID())
  {

    RealVector& linearData = _model->getPhysicalData();
    const CFreal gamma = _model->getgamma();
    const CFreal rho0  = linearData[LinEulerTerm::rho0];;
    const CFreal P0    = linearData[LinEulerTerm::P0];;
    const CFreal U0    = linearData[LinEulerTerm::U0];
    const CFreal V0    = linearData[LinEulerTerm::V0];

    const CFreal u = pdata[LinEulerTerm::u];
    const CFreal v = pdata[LinEulerTerm::v];
    const CFreal rho = pdata[LinEulerTerm::rho];
    const CFreal p = pdata[LinEulerTerm::p];

 EquationSetData& eqSet = LinEuler2DVarSet::getEqSetData()[0];
 const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();

 if (varIDs.size() == 4) {
      _physFlux(varIDs[0],XX) = U0*rho+u*rho0;
      _physFlux(varIDs[0],YY) = V0*rho+v*rho0;

      _physFlux(varIDs[1],XX) = rho0*U0*u+p;
      _physFlux(varIDs[1],YY) = rho0*V0*u;

      _physFlux(varIDs[2],XX) = rho0*U0*v;
      _physFlux(varIDs[2],YY) = rho0*V0*v+p;

      _physFlux(varIDs[3],XX) = U0*p+u*gamma*P0;
      _physFlux(varIDs[3],YY) = V0*p+v*gamma*P0;
    }
 else if (varIDs.size() == 3) {
      CF_DEBUG_EXIT;
      _physFlux(varIDs[0],XX) = U0*rho+u*rho0;
      _physFlux(varIDs[0],YY) = V0*rho+v*rho0;

      _physFlux(varIDs[1],XX) = rho0*U0*u+p;
      _physFlux(varIDs[1],YY) = rho0*V0*u;

      _physFlux(varIDs[2],XX) = rho0*U0*v;
      _physFlux(varIDs[2],YY) = rho0*V0*v+p;
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearizedEuler

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
