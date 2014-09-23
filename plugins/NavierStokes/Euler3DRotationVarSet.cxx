#include "NavierStokes/NavierStokes.hh"
#include "Euler3DRotationVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Euler3DRotationVarSet::Euler3DRotationVarSet(Common::SafePtr<BaseTerm> term) :
  EulerVarSet(term)
{
  // coordinates are needed in the physical data
  _physDataNeedCoordinates = true;
}
      
//////////////////////////////////////////////////////////////////////////////

Euler3DRotationVarSet::~Euler3DRotationVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRotationVarSet::setup()
{
  EulerVarSet::setup();

  // set EquationSetData
  // set EquationSetData
  Euler3DRotationVarSet::getEqSetData().resize(1);
  Euler3DRotationVarSet::getEqSetData()[0].setup(0,0,5);
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DRotationVarSet::getMaxEigenValue(const RealVector& data,
					       const RealVector& normal)
{ 
  /// @TODO check this in some paper ....
  const CFreal un = getNormalSpeed(data, normal);
  return un + data[EulerTerm::A];
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFreal Euler3DRotationVarSet::getMaxAbsEigenValue(const RealVector& data,
						  const RealVector& normal)
{
  // check this in some paper
  const CFreal un = getNormalSpeed(data, normal);
  return fabs(un) + data[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRotationVarSet::computeEigenValues(const RealVector& data,
					       const RealVector& normal,
					       RealVector& result)
{
  // check this in some paper
  const CFreal un = getNormalSpeed(data, normal);
  const vector<CFuint>& varIDs = Euler3DRotationVarSet::getEqSetData()[0].getEqSetVarIDs();
  
  result[varIDs[0]] = un;
  result[varIDs[1]] = un;
  result[varIDs[2]] = un;
  result[varIDs[3]] = un + data[EulerTerm::A];
  result[varIDs[4]] = un - data[EulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRotationVarSet::computeFlux (const RealVector& data,
					 const RealVector& normals)
{  
  const CFreal omega = getModel()->getOmega();
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
  const CFreal u = data[EulerTerm::VX];
  const CFreal v = data[EulerTerm::VY];
  const CFreal w = data[EulerTerm::VZ];
  const CFreal vn = u*nx + v*ny + w*nz;
  
  const CFreal perV =  - omega*data[EulerTerm::ZP];
  const CFreal perW =    omega*data[EulerTerm::YP];
  
  //  const CFreal relVn = vn - (perU*nx + perV*ny + perW*nz);  // 
  const CFreal relVn = vn - (perV*ny + perW*nz);  //  relative velocity * normal
  const CFreal rhoRelVn = data[EulerTerm::RHO]*relVn;
  const CFreal p = data[EulerTerm::P];
  
  const vector<CFuint>& varIDs = Euler3DRotationVarSet::getEqSetData()[0].getEqSetVarIDs();

  _fluxArray[varIDs[0]] = rhoRelVn;
  _fluxArray[varIDs[1]] = p*nx + u*rhoRelVn;
  _fluxArray[varIDs[2]] = p*ny + v*rhoRelVn;
  _fluxArray[varIDs[3]] = p*nz + w*rhoRelVn;
  _fluxArray[varIDs[4]] = rhoRelVn*data[EulerTerm::H] + p*(perV*ny + perW*nz);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRotationVarSet::computeStateFlux (const RealVector& data)
{
  // const RealVector& data = *_currData;
  
  //   const CFreal u = data[EulerTerm::VX];
  //   const CFreal v = data[EulerTerm::VY];
  //   const CFreal w = data[EulerTerm::VZ];
//   const CFreal rhoU = data[EulerTerm::RHO]*u;
//   const CFreal rhoV = data[EulerTerm::RHO]*v;
//   const CFreal rhoW = data[EulerTerm::RHO]*w;
//   const CFreal rhoUU = rhoU*u;
//   const CFreal rhoVV = rhoV*v;
//   const CFreal rhoWW = rhoW*w;
//   const CFreal rhoUV = rhoU*v;
//   const CFreal rhoUW = rhoU*w;
//   const CFreal rhoVW = rhoV*w;

//   _physFlux(0,XX) = rhoU;
//   _physFlux(0,YY) = rhoV;
//   _physFlux(0,ZZ) = rhoW;
  
//   _physFlux(1,XX) = data[EulerTerm::P] + rhoUU;
//   _physFlux(1,YY) = rhoUV;
//   _physFlux(1,ZZ) = rhoUW;
  
//   _physFlux(2,XX) = rhoUV;
//   _physFlux(2,YY) = data[EulerTerm::P] + rhoVV;
//   _physFlux(2,ZZ) = rhoVW;
  
//   _physFlux(3,XX) = rhoUW;
//   _physFlux(3,YY) = rhoVW;
//   _physFlux(3,ZZ) = data[EulerTerm::P] + rhoWW;
  
//   _physFlux(4,XX) = rhoU*data[EulerTerm::H];
//   _physFlux(4,YY) = rhoV*data[EulerTerm::H];
//   _physFlux(4,ZZ) = rhoW*data[EulerTerm::H];
  
  throw Common::NotImplementedException (FromHere(),"Euler3DRotationVarSet::computeFlux()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
