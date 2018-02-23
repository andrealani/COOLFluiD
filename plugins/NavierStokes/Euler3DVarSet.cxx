#include "NavierStokes/NavierStokes.hh"
#include "Euler3DVarSet.hh"
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

Euler3DVarSet::Euler3DVarSet(Common::SafePtr<BaseTerm> term) :
  EulerVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DVarSet::~Euler3DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DVarSet::setup()
{
  EulerVarSet::setup();

  // set EquationSetData
  Euler3DVarSet::getEqSetData().resize(1);
  Euler3DVarSet::getEqSetData()[0].setup(0,0,5);
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DVarSet::getMaxEigenValue(const RealVector& data,
				       const RealVector& normal)
{
  const CFreal un = getNormalSpeed(data, normal);
  return un + data[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DVarSet::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  const CFreal un = getNormalSpeed(data, normal);
  return fabs(un) + data[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DVarSet::computeEigenValues(const RealVector& data,
				       const RealVector& normal,
				       RealVector& result)
{
  const CFreal un = getNormalSpeed(data, normal);
  const vector<CFuint>& varIDs = Euler3DVarSet::getEqSetData()[0].getEqSetVarIDs();
  
  if (varIDs.size() == 5) {
    result[varIDs[0]] = un;
    result[varIDs[1]] = un;
    result[varIDs[2]] = un;
    result[varIDs[3]] = un + data[EulerTerm::A];
    result[varIDs[4]] = un - data[EulerTerm::A];
  }
  else if (varIDs.size() == 4) {
    result[varIDs[0]] = un;
    result[varIDs[1]] = un;
    result[varIDs[2]] = un + data[EulerTerm::A];
    result[varIDs[3]] = un - data[EulerTerm::A];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DVarSet::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = (normals.size() == 3) ? normals[ZZ] : 0.;
  const CFreal u = data[EulerTerm::VX];
  const CFreal v = data[EulerTerm::VY];
  const CFreal w = data[EulerTerm::VZ];
  const CFreal rhoVn = data[EulerTerm::RHO]*(u*nx + v*ny + w*nz);
  const vector<CFuint>& varIDs = Euler3DVarSet::getEqSetData()[0].getEqSetVarIDs();
  
  CFLog(DEBUG_MIN, "Euler2DVarSet::computeFlux() => " << 
	  CFPrintContainer<const vector<CFuint> >("varIDs = ", &varIDs));
  
  if (varIDs.size() == 5) {
    _fluxArray[varIDs[0]] = rhoVn;
    _fluxArray[varIDs[1]] = data[EulerTerm::P]*nx + u*rhoVn;
    _fluxArray[varIDs[2]] = data[EulerTerm::P]*ny + v*rhoVn;
    _fluxArray[varIDs[3]] = data[EulerTerm::P]*nz + w*rhoVn;
    _fluxArray[varIDs[4]] = rhoVn*data[EulerTerm::H];
  }
  else if (varIDs.size() == 4) {
    _fluxArray[varIDs[0]] = data[EulerTerm::P]*nx + u*rhoVn;
    _fluxArray[varIDs[1]] = data[EulerTerm::P]*ny + v*rhoVn;
    _fluxArray[varIDs[2]] = data[EulerTerm::P]*nz + w*rhoVn;
    _fluxArray[varIDs[3]] = rhoVn*data[EulerTerm::H];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DVarSet::computeStateFlux (const RealVector& data)
{
  // AL: this needs to be adapted for 2D and 1/2
  cf_assert(PhysicalModelStack::getActive()->getDim() == 3);
  
  const CFreal u = data[EulerTerm::VX];
  const CFreal v = data[EulerTerm::VY];
  const CFreal w = data[EulerTerm::VZ];
  const CFreal rhoU = data[EulerTerm::RHO]*u;
  const CFreal rhoV = data[EulerTerm::RHO]*v;
  const CFreal rhoW = data[EulerTerm::RHO]*w;
  const CFreal rhoUU = rhoU*u;
  const CFreal rhoVV = rhoV*v;
  const CFreal rhoWW = rhoW*w;
  const CFreal rhoUV = rhoU*v;
  const CFreal rhoUW = rhoU*w;
  const CFreal rhoVW = rhoV*w;
  const vector<CFuint>& varIDs = Euler3DVarSet::getEqSetData()[0].getEqSetVarIDs();

  if (varIDs.size() == 5) {
    _physFlux(0,XX) = rhoU;
    _physFlux(0,YY) = rhoV;
    _physFlux(0,ZZ) = rhoW;

    _physFlux(1,XX) = data[EulerTerm::P] + rhoUU;
    _physFlux(1,YY) = rhoUV;
    _physFlux(1,ZZ) = rhoUW;

    _physFlux(2,XX) = rhoUV;
    _physFlux(2,YY) = data[EulerTerm::P] + rhoVV;
    _physFlux(2,ZZ) = rhoVW;

    _physFlux(3,XX) = rhoUW;
    _physFlux(3,YY) = rhoVW;
    _physFlux(3,ZZ) = data[EulerTerm::P] + rhoWW;

    _physFlux(4,XX) = rhoU*data[EulerTerm::H];
    _physFlux(4,YY) = rhoV*data[EulerTerm::H];
    _physFlux(4,ZZ) = rhoW*data[EulerTerm::H];
  }
  else if (varIDs.size() == 4) {
    _physFlux(0,XX) = data[EulerTerm::P] + rhoUU;
    _physFlux(0,YY) = rhoUV;
    _physFlux(0,ZZ) = rhoUW;

    _physFlux(1,XX) = rhoUV;
    _physFlux(1,YY) = data[EulerTerm::P] + rhoVV;
    _physFlux(1,ZZ) = rhoVW;

    _physFlux(2,XX) = rhoUW;
    _physFlux(2,YY) = rhoVW;
    _physFlux(2,ZZ) = data[EulerTerm::P] + rhoWW;

    _physFlux(3,XX) = rhoU*data[EulerTerm::H];
    _physFlux(3,YY) = rhoV*data[EulerTerm::H];
    _physFlux(3,ZZ) = rhoW*data[EulerTerm::H];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
