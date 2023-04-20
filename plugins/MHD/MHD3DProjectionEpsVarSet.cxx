#include "MHD/MHD3DProjectionEpsVarSet.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsVarSet::MHD3DProjectionEpsVarSet(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DProjectionVarSet(term),
  _modelEps(term.d_castTo<MHDProjectionEpsTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsVarSet::~MHD3DProjectionEpsVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsVarSet::computeFlux (const RealVector& data,
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
    (rho*(u*u + v*v + w*w) + B2);
  
  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;
  
  _fluxArray[0] = Vn*rho;
  _fluxArray[1] = Vn*rho*u - Bn*Bx + P*nx;
  _fluxArray[2] = Vn*rho*v - Bn*By + P*ny;
  _fluxArray[3] = Vn*rho*w - Bn*Bz + P*nz;
  _fluxArray[4] = (v*Bx - By*u)*ny + (w*Bx - Bz*u)*nz + phi*nx;
  _fluxArray[5] = (u*By - Bx*v)*nx + (w*By - Bz*v)*nz + phi*ny;
  _fluxArray[6] = (u*Bz - Bx*w)*nx + (v*Bz - By*w)*ny + phi*nz;
  _fluxArray[7] = Vn*(E + P) - Bn*VdotB;
  _fluxArray[8] = refSpeedSq*Bn;
  
  // AL: new components to implement
  _fluxArray[9] = 0.;
  _fluxArray[10] = 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsVarSet::computeStateFlux (const RealVector& data)
{
  const CFreal rho = data[MHDTerm::RHO];
  const CFreal u = data[MHDTerm::VX];
  const CFreal v = data[MHDTerm::VY];
  const CFreal w = data[MHDTerm::VZ];
  const CFreal Bx = data[MHDTerm::BX];
  const CFreal By = data[MHDTerm::BY];
  const CFreal Bz = data[MHDTerm::BZ];
  const CFreal BxBy = Bx*By;
  const CFreal BxBz = Bx*Bz;
  const CFreal ByBz = By*Bz;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
  const CFreal VdotB = u*Bx  + v*By + w*Bz;
  const CFreal p = data[MHDTerm::P];
  const CFreal phi = data[MHDProjectionTerm::PHI];
  const CFreal P = p + 0.5*B2;
  const CFreal E = p/(getModel()->getGamma() - 1.0) + 0.5*
    (rho*(u*u + v*v + w*w) + B2);
  const CFreal rhoU = rho*u;
  const CFreal rhoV = rho*v;
  const CFreal rhoW = rho*w;
  const CFreal rhoUU = rhoU*u;
  const CFreal rhoVV = rhoV*v;
  const CFreal rhoWW = rhoW*w;
  const CFreal rhoUV = rhoU*v;
  const CFreal rhoUW = rhoU*w;
  const CFreal rhoVW = rhoV*w;
  const CFreal vBx = v*Bx;
  const CFreal wBx = w*Bx;
  const CFreal uBy = u*By;
  const CFreal wBy = w*By;
  const CFreal uBz = u*Bz;
  const CFreal vBz = v*Bz;

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _physFlux(0,XX) = rhoU;
  _physFlux(0,YY) = rhoV;
  _physFlux(0,ZZ) = rhoW;

  _physFlux(1,XX) = rhoUU + Bx*Bx + P;
  _physFlux(1,YY) = rhoUV + BxBy;
  _physFlux(1,ZZ) = rhoUW + BxBz;

  _physFlux(2,XX) = rhoUV + BxBy;
  _physFlux(2,YY) = rhoVV + By*By+ P;
  _physFlux(2,ZZ) = rhoVW + ByBz;

  _physFlux(3,XX) = rhoUW + BxBz;
  _physFlux(3,YY) = rhoVW + ByBz;
  _physFlux(3,ZZ) = rhoWW + Bz*Bz + P;

  _physFlux(4,XX) = phi;
  _physFlux(4,YY) = vBx - uBy;
  _physFlux(4,ZZ) = wBx - uBz;

  _physFlux(5,XX) = uBy - vBx;
  _physFlux(5,YY) = phi;
  _physFlux(5,ZZ) = wBy - vBz;

  _physFlux(6,XX) = uBz - wBx;
  _physFlux(6,YY) = vBz - wBy;
  _physFlux(6,ZZ) = phi;

  _physFlux(7,XX) = u*(E + P) - Bx*VdotB;
  _physFlux(7,YY) = v*(E + P) - By*VdotB;
  _physFlux(7,ZZ) = w*(E + P) - Bz*VdotB;

  _physFlux(8,XX) = refSpeedSq*Bx;
  _physFlux(8,YY) = refSpeedSq*By;
  _physFlux(8,ZZ) = refSpeedSq*Bz;

  // AL: new components to implement
  _physFlux(9,XX) = 0.;
  _physFlux(9,YY) = 0.;
  _physFlux(9,ZZ) = 0.;
  
  _physFlux(10,XX) = 0.;
  _physFlux(10,YY) = 0.;
  _physFlux(10,ZZ) = 0.;
}

//////////////////////////////////////////////////////////////////////
 
CFreal MHD3DProjectionEpsVarSet::getMaxEigenValue(const RealVector& data,
						  const RealVector& normal)
{
  // AL: included B0 contribution which was missing
  CFreal Bx0 = 0.0, By0 = 0.0, Bz0 = 0.0;
  const CFreal Bx1 = data[MHDTerm::BX];
  const CFreal By1 = data[MHDTerm::BY];
  const CFreal Bz1 = data[MHDTerm::BZ];
  const CFreal BxTotal = Bx1 + Bx0;
  const CFreal ByTotal = By1 + By0;
  const CFreal BzTotal = Bz1 + Bz0;
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + BzTotal*BzTotal;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY] + data[MHDTerm::VZ]*normal[ZZ];
  const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY] + BzTotal*normal[ZZ];
  const CFreal sqrbar = sqrt(data[MHDTerm::RHO]);
  const CFreal astar2 = (getModel()->getGamma()*p + sqBTotal)*invRho;
  const CFreal astarb = sqrt(astar2*astar2 - 4.0*getModel()->getGamma()*p*BnTotal*BnTotal*invRho*invRho);
  const CFreal cf2 = 0.5*(astar2 + astarb);
  // const CFreal ca = std::abs(BnTotal)/sqrbar;
  const CFreal cf = sqrt(cf2);
  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal maxEigenValue = max(refSpeed,(Vn + cf));
  return maxEigenValue;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DProjectionEpsVarSet::getMaxAbsEigenValue(const RealVector& data,
						  const RealVector& normal)
{
   // AL: included B0 contribution which was missing
  const std::string& potentialBType = getModel()->getPotentialBType();
  
  CFreal Bx0 = 0.0, By0 = 0.0, Bz0 = 0.0;
  const CFreal BxTotal = data[MHDTerm::BX] + Bx0;
  const CFreal ByTotal = data[MHDTerm::BY] + By0;
  const CFreal BzTotal = data[MHDTerm::BZ] + Bz0;
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + BzTotal*BzTotal;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY] + data[MHDTerm::VZ]*normal[ZZ];
  const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY] + BzTotal*normal[ZZ];
  const CFreal sqrbar = sqrt(data[MHDTerm::RHO]);
  const CFreal astar2 = (getModel()->getGamma()*p + sqBTotal)*invRho;
  const CFreal astarb = sqrt(astar2*astar2 - 4.0*getModel()->getGamma()*p*BnTotal*BnTotal*invRho*invRho);
  const CFreal cf2 = 0.5*(astar2 + astarb);
  // const CFreal ca = std::abs(BnTotal)/sqrbar;
  const CFreal cf = sqrt(cf2);
  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal maxEigenValue = max(refSpeed,(std::abs(Vn) + cf));
  return maxEigenValue;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsVarSet::computeEigenValues(const RealVector& data,
					       const RealVector& normal,
					       RealVector& result)
{
  // eigenvalues for the flows not involving any background potential magnetic field are also correctly
  // calculated in this function since B0 field will automatically be zero in that case
  // eigenvalues based on total magnetic field, Btotal
  
  CFreal Bx0 = 0.0, By0 = 0.0, Bz0 = 0.0;
  const CFreal BxTotal = data[MHDTerm::BX] + Bx0;
  const CFreal ByTotal = data[MHDTerm::BY] + By0;
  const CFreal BzTotal = data[MHDTerm::BZ] + Bz0;
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + BzTotal*BzTotal;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY] + data[MHDTerm::VZ]*normal[ZZ];
  const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY] + BzTotal*normal[ZZ];
  const CFreal sqrbar = sqrt(data[MHDTerm::RHO]);
  const CFreal astar2 = (getModel()->getGamma()*p + sqBTotal)*invRho;
  const CFreal astarb = sqrt(astar2*astar2 - 4.0*getModel()->getGamma()*p*BnTotal*BnTotal*invRho*invRho);
  const CFreal cf2 = 0.5*(astar2 + astarb);
  const CFreal cs2 = 0.5*(astar2 - astarb);
  const CFreal ca = std::abs(BnTotal)/sqrbar;
  const CFreal cf = sqrt(cf2);
  const CFreal cs = sqrt(cs2);
  const CFreal refSpeed = getModel()->getRefSpeed();
  
  result[0] = Vn - cf;
  result[1] = Vn - ca;
  result[2] = Vn - cs;
  result[3] = Vn;
  result[4] = refSpeed;
  result[5] = -refSpeed;
  result[6] = Vn + cs;
  result[7] = Vn + ca;
  result[8] = Vn + cf;

  // AL: new components to implement
  result[9] = 0.;
  result[10] = 0.;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
