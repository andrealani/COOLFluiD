#include "NavierStokes/NavierStokes.hh"
#include "Euler2DVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFPrintContainer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Euler2DVarSet::Euler2DVarSet(Common::SafePtr<BaseTerm> term) :
  EulerVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DVarSet::~Euler2DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DVarSet::setup()
{
  EulerVarSet::setup();

  // set EquationSetData
  Euler2DVarSet::getEqSetData().resize(1);
  Euler2DVarSet::getEqSetData()[0].setup(0,0,4);
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DVarSet::getMaxEigenValue(const RealVector& pdata,
				       const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return un + pdata[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
					  const RealVector& normal)
{
  const CFreal un = getNormalSpeed(pdata, normal);
  return fabs(un) + pdata[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DVarSet::computeEigenValues(const RealVector& pdata,
				       const RealVector& normal,
				       RealVector& result)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  
  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID()) {
    const vector<CFuint>& varIDs = Euler2DVarSet::getEqSetData()[0].getEqSetVarIDs();
    const CFreal un = getNormalSpeed(pdata, normal);
    if (varIDs.size() == 4) {
      
      result[varIDs[0]] = un;
      result[varIDs[1]] = un;
      result[varIDs[2]] = un + pdata[EulerTerm::A];
      result[varIDs[3]] = un - pdata[EulerTerm::A];
    }
    else if (varIDs.size() == 3) {
      
      result[varIDs[0]] = un;
      result[varIDs[1]] = un + pdata[EulerTerm::A];
      result[varIDs[2]] = un - pdata[EulerTerm::A];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DVarSet::computeFlux (const RealVector& pdata, const RealVector& normals)
{
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
   getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID()) {
    const CFreal nx = normals[XX];
    const CFreal ny = normals[YY];
    const CFreal u = pdata[EulerTerm::VX];
    const CFreal v = pdata[EulerTerm::VY];
    const CFreal un = getNormalSpeed(pdata, normals);
    const CFreal rhoVn = pdata[EulerTerm::RHO]*un;
    
    EquationSetData& eqSet = Euler2DVarSet::getEqSetData()[0];
    const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();
    
    CFLog(DEBUG_MIN, "Euler2DVarSet::computeFlux() => " << 
	  CFPrintContainer<const vector<CFuint> >("varIDs = ", &varIDs));
    
    if (varIDs.size() == 4) {
      _fluxArray[varIDs[0]] = rhoVn;
      _fluxArray[varIDs[1]] = pdata[EulerTerm::P]*nx + u*rhoVn;
      _fluxArray[varIDs[2]] = pdata[EulerTerm::P]*ny + v*rhoVn;
      _fluxArray[varIDs[3]] = rhoVn*pdata[EulerTerm::H];
    }
    else if (varIDs.size() == 3) {
      _fluxArray[varIDs[0]] = pdata[EulerTerm::P]*nx + u*rhoVn;
      _fluxArray[varIDs[1]] = pdata[EulerTerm::P]*ny + v*rhoVn;
      _fluxArray[varIDs[2]] = rhoVn*pdata[EulerTerm::H];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DVarSet::computeStateFlux (const RealVector& pdata)
{  
  EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->
    getEquationSubSysDescriptor();
  const CFuint totalNbEqSS = eqSS.getNbEqsSS();
  if (totalNbEqSS == _nbEqs || eqSS.getEqSS() == getEqSetData()[0].getEqSetID()) {
    const CFreal u = pdata[EulerTerm::VX];
    const CFreal v = pdata[EulerTerm::VY];
    const CFreal rhoU = pdata[EulerTerm::RHO]*u;
    const CFreal rhoV = pdata[EulerTerm::RHO]*v;
    const CFreal rhoUU = rhoU*u;
    const CFreal rhoUV = rhoU*v;
    const CFreal rhoVV = rhoV*v;
    
    EquationSetData& eqSet = Euler2DVarSet::getEqSetData()[0];
    const vector<CFuint>& varIDs = eqSet.getEqSetVarIDs();
    
    if (varIDs.size() == 4) {
      _physFlux(varIDs[0],XX) = rhoU;
      _physFlux(varIDs[0],YY) = rhoV;
      
      _physFlux(varIDs[1],XX) = pdata[EulerTerm::P] + rhoUU;
      _physFlux(varIDs[1],YY) = rhoUV;
      
      _physFlux(varIDs[2],XX) = rhoUV;
      _physFlux(varIDs[2],YY) = pdata[EulerTerm::P] + rhoVV;
      
      _physFlux(varIDs[3],XX) = rhoU*pdata[EulerTerm::H];
      _physFlux(varIDs[3],YY) = rhoV*pdata[EulerTerm::H];
    }
    else if (varIDs.size() == 3) {
      _physFlux(varIDs[0],XX) = pdata[EulerTerm::P] + rhoUU;
      _physFlux(varIDs[0],YY) = rhoUV;
      
      _physFlux(varIDs[1],XX) = rhoUV;
      _physFlux(varIDs[1],YY) = pdata[EulerTerm::P] + rhoVV;
      
      _physFlux(varIDs[2],XX) = rhoU*pdata[EulerTerm::H];
      _physFlux(varIDs[2],YY) = rhoV*pdata[EulerTerm::H];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

#if 0 // keep this function here untill we implement the FluctSplit MetaSchemes properly
void Euler2DVarSet::computeFlux ( const State& vars, const RealVector& physdata, RealMatrix& flux) const
{
  const CFreal u     = physdata[EulerTerm::VX];
  const CFreal v     = physdata[EulerTerm::VY];
  const CFreal rhoU  = physdata[EulerTerm::RHO]*u;
  const CFreal rhoV  = physdata[EulerTerm::RHO]*v;
  const CFreal rhoUU = rhoU*u;
  const CFreal rhoUV = rhoU*v;
  const CFreal rhoVV = rhoV*v;

  flux(0,XX) = rhoU;
  flux(0,YY) = rhoV;

  flux(1,XX) = physdata[EulerTerm::P] + rhoUU;
  flux(1,YY) = rhoUV;

  flux(2,XX) = rhoUV;
  flux(2,YY) = physdata[EulerTerm::P] + rhoVV;

  flux(3,XX) = rhoU*physdata[EulerTerm::H];
  flux(3,YY) = rhoV*physdata[EulerTerm::H];
}
#endif

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
