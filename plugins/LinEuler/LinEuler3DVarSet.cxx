#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler3DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

LinEuler3DVarSet::LinEuler3DVarSet(Common::SafePtr<BaseTerm> term) :
  LinEulerVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

LinEuler3DVarSet::~LinEuler3DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DVarSet::setup()
{
  LinEulerVarSet::setup();

  // set EquationSetData
  LinEuler3DVarSet::getEqSetData().resize(1);
  LinEuler3DVarSet::getEqSetData()[0].setup(0,0,5);
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinEuler3DVarSet::getMaxEigenValue(const RealVector& pdata,
               const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return un + pdata[LinEulerTerm::c];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LinEuler3DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
            const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return fabs(un) + pdata[LinEulerTerm::c];
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DVarSet::computeEigenValues(const RealVector& pdata,
               const RealVector& normal, RealVector& result)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();

  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID())
  {

   const vector<CFuint>& varIDs = LinEuler3DVarSet::getEqSetData()[0].getEqSetVarIDs();
   const CFreal un = getNormalSpeed(pdata, normal);

      result[varIDs[0]] = un;
      result[varIDs[1]] = un;
      result[varIDs[2]] = un;
      result[varIDs[3]] = un + pdata[LinEulerTerm::c];
      result[varIDs[4]] = un - pdata[LinEulerTerm::c];

  }
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DVarSet::computeFlux (const RealVector& pdata, const RealVector& normals)
{

  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();

  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID())
  {
    RealVector& linearData = _model->getPhysicalData();
    const CFreal gamma = linearData[LinEulerTerm::GAMMA];
    const CFreal rho0  = linearData[LinEulerTerm::rho0];;
    const CFreal P0    = linearData[LinEulerTerm::P0];;
    const CFreal U0    = linearData[LinEulerTerm::U0];
    const CFreal V0    = linearData[LinEulerTerm::V0];
    const CFreal W0    = linearData[LinEulerTerm::W0];

    const CFreal nx = normals[XX];
    const CFreal ny = normals[YY];
    const CFreal nz = normals[ZZ];
    const CFreal Vn = U0*nx+V0*ny+W0*nz;
    const CFreal u = pdata[LinEulerTerm::u];
    const CFreal v = pdata[LinEulerTerm::v];
    const CFreal w = pdata[LinEulerTerm::w];
    const CFreal un = u*nx + v*ny + w*nz;
    const CFreal rho = pdata[LinEulerTerm::rho];
    const CFreal p = pdata[LinEulerTerm::p];

 EquationSetData& eqSet = LinEuler3DVarSet::getEqSetData()[0];
 const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();

     _fluxArray[varIDs[0]] = Vn*rho+un*rho0;
     _fluxArray[varIDs[1]] = Vn*rho0*u+p*nx;
     _fluxArray[varIDs[2]] = Vn*rho0*v+p*ny;
     _fluxArray[varIDs[3]] = Vn*rho0*w+p*nz;
     _fluxArray[varIDs[4]] = Vn*p+un*gamma*P0;

  }
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DVarSet::computeStateFlux (const RealVector& pdata)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID())
  {
    RealVector& linearData = _model->getPhysicalData();
    const CFreal gamma = linearData[LinEulerTerm::GAMMA];
    const CFreal rho0  = linearData[LinEulerTerm::rho0 ];
    const CFreal P0    = linearData[LinEulerTerm::P0   ];
    const CFreal U0    = linearData[LinEulerTerm::U0   ];
    const CFreal V0    = linearData[LinEulerTerm::V0   ];
    const CFreal W0    = linearData[LinEulerTerm::W0   ];

    const CFreal u = pdata[LinEulerTerm::u];
    const CFreal v = pdata[LinEulerTerm::v];
    const CFreal w = pdata[LinEulerTerm::w];
    const CFreal rho = pdata[LinEulerTerm::rho];
    const CFreal p = pdata[LinEulerTerm::p];

    EquationSetData& eqSet = LinEuler3DVarSet::getEqSetData()[0];
    const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();

    if (varIDs.size() == 5) {
      _physFlux(varIDs[0],XX) = U0*rho + u*rho0;
      _physFlux(varIDs[0],YY) = V0*rho + v*rho0;
      _physFlux(varIDs[0],ZZ) = W0*rho + w*rho0;

      _physFlux(varIDs[1],XX) = rho0*U0*u + p;
      _physFlux(varIDs[1],YY) = rho0*V0*u;
      _physFlux(varIDs[1],ZZ) = rho0*W0*u;

      _physFlux(varIDs[2],XX) = rho0*U0*v;
      _physFlux(varIDs[2],YY) = rho0*V0*v + p;
      _physFlux(varIDs[2],ZZ) = rho0*W0*v;

      _physFlux(varIDs[3],XX) = rho0*U0*w;
      _physFlux(varIDs[3],YY) = rho0*V0*w;
      _physFlux(varIDs[3],ZZ) = rho0*W0*w + p;

      _physFlux(varIDs[4],XX) = U0*p + u*gamma*P0;
      _physFlux(varIDs[4],YY) = V0*p + v*gamma*P0;
      _physFlux(varIDs[4],ZZ) = W0*p + w*gamma*P0;
    }
    else if (varIDs.size() == 4) {
      CF_DEBUG_EXIT;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearizedEuler

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
