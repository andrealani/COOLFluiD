#include "MHD/MHD.hh"

#include "MHD/MHD3DProjectionPrimE.hh"
#include "Common/CFLog.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionPrimE, ConvectiveVarSet, MHDModule, 1> 
mhd3DProjectionPrimEProvider("MHD3DProjectionPrimE");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPrimE::MHD3DProjectionPrimE(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DProjectionPrim(term)
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPrimE::~MHD3DProjectionPrimE()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::setup()
{
  MHD3DProjectionVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD3DProjectionPrimE::getExtraVarNames() const
{
  return MHD3DProjectionPrim::getExtraVarNames();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionPrimE::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::computeEigenValuesVectors(RealMatrix& rightEv,
                                      RealMatrix& leftEv,
                                      RealVector& eValues,
                                      const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionPrimE::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::setDimensionalValuesPlusExtraValues(const State& state,
                                                              RealVector& result,
                                                              RealVector& extra)
{
  MHD3DProjectionPrim::setDimensionalValuesPlusExtraValues(state, result, extra); 
  
  // override the total rhoE entry
  const CFreal gammaMinus1 = getModel()->getGamma() - 1.;
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal V2 = u*u + v*v + w*w;
  extra[7] = (state[7]/gammaMinus1) + 0.5*(rho*V2);
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::splitJacobian(RealMatrix& jacobPlus,
					 RealMatrix& jacobMin,
					 RealVector& eValues,
					 const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionPrimE::splitJacobian()"); 
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::computePhysicalData(const State& state,
					      RealVector& data)
{
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;

  data[MHDProjectionTerm::VX] = u;
  data[MHDProjectionTerm::VY] = v;
  data[MHDProjectionTerm::VZ] = w;
  data[MHDProjectionTerm::BX] = Bx;
  data[MHDProjectionTerm::BY] = By;
  data[MHDProjectionTerm::BZ] = Bz;
  data[MHDProjectionTerm::RHO] = rho;
  data[MHDProjectionTerm::V] = sqrt(V2);
  data[MHDProjectionTerm::B] = sqrt(B2);
  data[MHDProjectionTerm::P] = state[7];
  data[MHDProjectionTerm::PHI] = state[8];

  // Haopeng: modify speed of sound "A" if needed
  data[MHDProjectionTerm::A] = sqrt(getModel()->getGamma()*data[MHDProjectionTerm::P]/rho);
  
  data[MHDProjectionTerm::GAMMA] = getModel()->getGamma();
  
  const RealVector& node = state.getCoordinates();
  data[MHDTerm::XP] =  node[XX];
  data[MHDTerm::YP] =  node[YY];
  data[MHDTerm::ZP] =  node[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::computeStateFromPhysicalData(const RealVector& data,
						  State& state)
{
  state[0] = data[MHDProjectionTerm::RHO];
  state[1] = data[MHDProjectionTerm::VX];
  state[2] = data[MHDProjectionTerm::VY];
  state[3] = data[MHDProjectionTerm::VZ];
  state[4] = data[MHDProjectionTerm::BX];
  state[5] = data[MHDProjectionTerm::BY];
  state[6] = data[MHDProjectionTerm::BZ];
  state[7] = data[MHDProjectionTerm::P];
  state[8] = data[MHDProjectionTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::computeFlux (const RealVector& data,
					const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
  const CFreal rho = data[MHDTerm::RHO];
  const CFreal u = data[MHDTerm::VX];
  const CFreal v = data[MHDTerm::VY];
  const CFreal w = data[MHDTerm::VZ];
  const CFreal Bx = data[MHDTerm::BX];
  const CFreal By = data[MHDTerm::BY];
  const CFreal Bz = data[MHDTerm::BZ];
  const CFreal Vn = u*nx  + v*ny + w*nz;
  const CFreal Bn = Bx*nx + By*ny + Bz*nz;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
  const CFreal VdotB = u*Bx  + v*By + w*Bz;
  const CFreal p = data[MHDTerm::P];
  const CFreal phi = data[MHDProjectionTerm::PHI];
  const CFreal P = p + 0.5*B2;
  const CFreal E = p/(getModel()->getGamma() - 1.0) + 0.5*
    (rho*(u*u + v*v + w*w));

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  // Haopeng: check this
  _fluxArray[0] = Vn*rho;
  _fluxArray[1] = Vn*rho*u - Bn*Bx + P*nx;
  _fluxArray[2] = Vn*rho*v - Bn*By + P*ny;
  _fluxArray[3] = Vn*rho*w - Bn*Bz + P*nz;
  _fluxArray[4] = (v*Bx - By*u)*ny + (w*Bx - Bz*u)*nz + phi*nx;
  _fluxArray[5] = (u*By - Bx*v)*nx + (w*By - Bz*v)*nz + phi*ny;
  _fluxArray[6] = (u*Bz - Bx*w)*nx + (v*Bz - By*w)*ny + phi*nz;
  _fluxArray[7] = Vn*(E + P);
  _fluxArray[8] = refSpeedSq*Bn;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrimE::computeStateFlux(const RealVector& data)
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionPrimE::computeStateFlux()"); 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

