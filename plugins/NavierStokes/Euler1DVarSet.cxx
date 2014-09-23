#include "NavierStokes/NavierStokes.hh"
#include "Euler1DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Euler1DVarSet::Euler1DVarSet(Common::SafePtr<BaseTerm> term) :
  EulerVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler1DVarSet::~Euler1DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DVarSet::setup()
{
  EulerVarSet::setup();

  // set EquationSetData
  Euler1DVarSet::getEqSetData().resize(1);
  Euler1DVarSet::getEqSetData()[0].setup(0,0,3);
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFreal Euler1DVarSet::getMaxEigenValue(const RealVector& pdata, 
				       const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return un + pdata[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler1DVarSet::getMaxAbsEigenValue(const RealVector& pdata, 
					  const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return fabs(un) + pdata[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DVarSet::computeEigenValues(const RealVector& pdata, 
				       const RealVector& normal, 
				       RealVector& eValues)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  
  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID()) {
    const vector<CFuint>& varIDs = Euler1DVarSet::getEqSetData()[0].getEqSetVarIDs();
    const CFreal un = getNormalSpeed(pdata, normal);;
    if (varIDs.size() == 3) {
      eValues[varIDs[0]] = un;
      eValues[varIDs[1]] = un + pdata[EulerTerm::A];
      eValues[varIDs[2]] = un - pdata[EulerTerm::A];
    }
    else if (varIDs.size() == 2) {
      eValues[varIDs[0]] = un + pdata[EulerTerm::A];
      eValues[varIDs[1]] = un - pdata[EulerTerm::A];
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DVarSet::computeFlux (const RealVector& pdata, const RealVector& normals)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  
  if (totalNbEqSS ==_nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID()) {
    const CFreal un = getNormalSpeed(pdata, normals);
    const CFreal nx = normals[XX];
    const CFreal rhoVn = pdata[EulerTerm::RHO]*un;
    const CFreal u = pdata[EulerTerm::VX];
    
    EquationSetData& eqSet = Euler1DVarSet::getEqSetData()[0];
    const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();

    if (varIDs.size() == 3) {
      _fluxArray[varIDs[0]] = rhoVn;
      _fluxArray[varIDs[1]] = pdata[EulerTerm::P]*nx + u*rhoVn;
      _fluxArray[varIDs[2]] = rhoVn*pdata[EulerTerm::H];
    }
    else if (varIDs.size() == 2) {
      _fluxArray[varIDs[0]] = pdata[EulerTerm::P]*nx + u*rhoVn;
      _fluxArray[varIDs[1]] = rhoVn*pdata[EulerTerm::H];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DVarSet::computeStateFlux (const RealVector& pdata)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  
  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID()) {
    const CFreal u = pdata[EulerTerm::VX];
    const CFreal rhoU = pdata[EulerTerm::RHO]*u;
    const CFreal rhoUU = rhoU*u;
    
    EquationSetData& eqSet = Euler1DVarSet::getEqSetData()[0];
    const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();
    
    if (varIDs.size() == 3) {
      _physFlux(varIDs[0],XX) = rhoU;
      _physFlux(varIDs[1],XX) = pdata[EulerTerm::P] + rhoUU;
      _physFlux(varIDs[2],XX) = rhoU*pdata[EulerTerm::H];
    }
    else if (varIDs.size() == 2) {
      _physFlux(varIDs[0],XX) = pdata[EulerTerm::P] + rhoUU;
      _physFlux(varIDs[1],XX) = rhoU*pdata[EulerTerm::H];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
