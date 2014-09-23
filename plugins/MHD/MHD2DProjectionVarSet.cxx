#include "MHD2DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionVarSet::MHD2DProjectionVarSet(Common::SafePtr<BaseTerm> term) :
  ConvectiveVarSet(term),
  _model(term.d_castTo<MHDProjectionTerm>()),
  _BDipole(),
  _tanakaFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionVarSet::~MHD2DProjectionVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionVarSet::computeFlux (const RealVector& data, 
					 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal rho = data[MHDTerm::RHO];
  const CFreal u = data[MHDTerm::VX];
  const CFreal v = data[MHDTerm::VY];
  const CFreal Bx = data[MHDTerm::BX];
  const CFreal By = data[MHDTerm::BY];
  const CFreal Bz = data[MHDTerm::BZ];
  const CFreal Vn = u*nx  + v*ny;
  const CFreal Bn = Bx*nx + By*ny;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
  const CFreal VdotB = u*Bx  + v*By;
  const CFreal p = data[MHDTerm::P];
  const CFreal phi = data[MHDProjectionTerm::PHI];
  const CFreal P = p + 0.5*B2;
  const CFreal E = p/(getModel()->getGamma() - 1.0) + 0.5*
    (rho*(u*u + v*v) + B2);

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _fluxArray[0] = Vn*rho;
  _fluxArray[1] = Vn*rho*u - Bn*Bx + P*nx;
  _fluxArray[2] = Vn*rho*v - Bn*By + P*ny;
  _fluxArray[3] = - Bn*Bz;
  _fluxArray[4] = (v*Bx - By*u)*ny + phi*nx;
  _fluxArray[5] = (u*By - Bx*v)*nx + phi*ny;
  _fluxArray[6] = u*Bz*nx + v*Bz*ny;
  _fluxArray[7] = Vn*(E + P) - Bn*VdotB;
  _fluxArray[8] = refSpeedSq*Bn;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHD2DProjectionVarSet::computeTanakaFluxPowell99Formulation(const RealVector& data,
                                                                        const RealVector& normals)
{
  // According to the explanations in Powell, K.G. et. al., JCP, Vol.154, pp.284-309, 1999

  computeMagneticDipole(data[MHDTerm::XP],data[MHDTerm::YP]);
  const CFreal Bx0 = _BDipole[0];
  const CFreal By0 = _BDipole[1];

  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal rho = data[MHDTerm::RHO];
  const CFreal u = data[MHDTerm::VX];
  const CFreal v = data[MHDTerm::VY];
  const CFreal w = data[MHDTerm::VZ];
  const CFreal Bx1 = data[MHDTerm::BX];
  const CFreal By1 = data[MHDTerm::BY];
  const CFreal Bz1 = data[MHDTerm::BZ];

  // n stands for normal to the face of the control volume

  const CFreal Vn = u*nx + v*ny;
  const CFreal Bn0 = Bx0*nx + By0*ny;

  const CFreal Bn1 = Bx1*nx + By1*ny;
  const CFreal sqB1 = Bx1*Bx1 + By1*By1 + Bz1*Bz1;
  const CFreal VdotB1 = u*Bx1  + v*By1 + w*Bz1;
  const CFreal B1dotB0 = Bx1*Bx0 + By1*By0;
  const CFreal p = data[MHDTerm::P];

  const CFreal P = p + 0.5*sqB1 + B1dotB0;
  const CFreal E1 = p/(getModel()->getGamma() - 1.)
          + 0.5 * (rho*(u*u + v*v + w*w) + sqB1);

  const CFreal phi = data[MHDProjectionTerm::PHI];

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _tanakaFlux[0] = Vn*rho;
  _tanakaFlux[1] = Vn*rho*u - Bn1*Bx1 - Bn0*Bx1 - Bn1*Bx0 + P*nx;
  _tanakaFlux[2] = Vn*rho*v - Bn1*By1 - Bn0*By1 - Bn1*By0 + P*ny;
  _tanakaFlux[3] = Vn*rho*w - Bn1*Bz1 - Bn0*Bz1;
  _tanakaFlux[4] = (v*Bx1 - By1*u)*ny + (v*Bx0 - By0*u)*ny + phi*nx;
  _tanakaFlux[5] = (u*By1 - Bx1*v)*nx + (u*By0 - Bx0*v)*nx + phi*ny;
  _tanakaFlux[6] = (u*Bz1 - Bx1*w)*nx + (v*Bz1 - By1*w)*ny + (Bx0*w)*nx + (By0*w)*ny;
  _tanakaFlux[7] = Vn*(E1 + p + 0.5*sqB1) - Bn1*VdotB1 + Vn*B1dotB0 - Bn0*VdotB1;
  _tanakaFlux[8] = refSpeedSq*Bn1;

  return _tanakaFlux;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionVarSet::computeStateFlux (const RealVector& data)
{
  const CFreal rho = data[MHDTerm::RHO];
  const CFreal u = data[MHDTerm::VX];
  const CFreal v = data[MHDTerm::VY];
  const CFreal Bx = data[MHDTerm::BX];
  const CFreal By = data[MHDTerm::BY];
  const CFreal Bz = data[MHDTerm::BZ];
  const CFreal BxBy = Bx*By;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
  const CFreal VdotB = u*Bx  + v*By;
  const CFreal vBx = v*Bx;
  const CFreal uBy = u*By;
  const CFreal p = data[MHDTerm::P];
  const CFreal phi = data[MHDProjectionTerm::PHI];
  const CFreal P = p + 0.5*B2;
  const CFreal E = p/(getModel()->getGamma() - 1.0) + 0.5*
    (rho*(u*u + v*v) + B2);
  const CFreal rhoU = rho*u;
  const CFreal rhoV = rho*v;
  const CFreal rhoUV = rhoU*v;
  const CFreal rhoUU = rhoU*u;
  const CFreal rhoVV = rhoV*v;

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _physFlux(0,XX) = rhoU;
  _physFlux(0,YY) = rhoV;

  _physFlux(1,XX) = rhoUU + Bx*Bx + P;
  _physFlux(1,YY) = rhoUV + BxBy;

  _physFlux(2,XX) = rhoUV + BxBy;
  _physFlux(2,YY) = rhoVV + By*By + P;

  _physFlux(3,XX) = -Bz*Bx;
  _physFlux(3,YY) = -Bz*By;

  _physFlux(4,XX) = phi;
  _physFlux(4,YY) = vBx - uBy;

  _physFlux(5,XX) = uBy - vBx;
  _physFlux(5,YY) = phi;

  _physFlux(6,XX) = u*Bz;
  _physFlux(6,YY) = v*Bz;

  _physFlux(7,XX) = u*(E + P) - Bx*VdotB;
  _physFlux(7,YY) = v*(E + P) - By*VdotB;

  _physFlux(8,XX) = refSpeedSq*Bx;
  _physFlux(8,YY) = refSpeedSq*By;
}

//////////////////////////////////////////////////////////////////////

CFreal MHD2DProjectionVarSet::getMaxEigenValue(const RealVector& data,
					       const RealVector& normal)
{
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal B = data[MHDTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY];
  const CFreal Bn = data[MHDTerm::BX]*normal[XX] + data[MHDTerm::BY]*normal[YY];
  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
				  - 4.0*getModel()->getGamma()*p*Bn*Bn*
				  invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal maxEigenValue = max(refSpeed,(Vn + cf));
  return maxEigenValue;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD2DProjectionVarSet::getMaxAbsEigenValue(const RealVector& data,
						  const RealVector& normal)
{
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal B = data[MHDTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY];
  const CFreal Bn = data[MHDTerm::BX]*normal[XX] + data[MHDTerm::BY]*normal[YY];
  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
          - 4.0*getModel()->getGamma()*p*Bn*Bn*
          invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal maxEigenValue = max(refSpeed,(std::abs(Vn) + cf));
  return maxEigenValue;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionVarSet::computeMagneticDipole(CFreal xCoord, CFreal yCoord)
{
  // this formula is valid for a dipole at the origin
  const CFreal mX = getMX();
  const CFreal mY = getMY();
  
  const CFreal r = sqrt(xCoord*xCoord + yCoord*yCoord);
  const CFreal rCubeInv = 1.0/(r*r*r);
  const CFreal rSqInv = 1.0/(r*r);
  const CFreal mDotr = mX*xCoord + mY*yCoord;

  const CFreal BxDipole = rCubeInv*(3.0*xCoord*rSqInv*mDotr - mX);
  const CFreal ByDipole = rCubeInv*(3.0*yCoord*rSqInv*mDotr - mY);

  _BDipole[0] = BxDipole;
  _BDipole[1] = ByDipole;
}

//////////////////////////////////////////////////////////////////////

void MHD2DProjectionVarSet::computeEigenValues(const RealVector& data,
					       const RealVector& normal,
					       RealVector& result)
{
  // eigenvalues for the flows not involving any dipole field are also correctly
  // calculated in this function since B0 field will automatically be zero in that case
  // eigenvalues based on total magnetic field, Btotal
  computeMagneticDipole(data[MHDTerm::XP], data[MHDTerm::YP]);
  const CFreal Bx0 = _BDipole[0];
  const CFreal By0 = _BDipole[1];
  
  const CFreal Bx1 = data[MHDTerm::BX];
  const CFreal By1 = data[MHDTerm::BY];
  const CFreal BxTotal = Bx1 + Bx0;
  const CFreal ByTotal = By1 + By0;
  const CFreal Bz1 = data[MHDTerm::BZ];

  const CFreal invRho = 1./data[MHDTerm::RHO];
  //const CFreal B = data[MHDTerm::B];
  const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + Bz1*Bz1;
  //const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] +
    data[MHDTerm::VY]*normal[YY];
  // const CFreal Bn = data[MHDTerm::BX]*normal[XX] +
//     data[MHDTerm::BY]*normal[YY];
  const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY];
  const CFreal sqrbar = sqrt(data[MHDTerm::RHO]);
  const CFreal astar2 = (getModel()->getGamma()*p + sqBTotal)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
				  - 4.0*getModel()->getGamma()*p*BnTotal*BnTotal*invRho*invRho));
  CFreal cs2 = 0.5*(astar2 - sqrt(astar2*astar2
				  - 4.0*getModel()->getGamma()*p*BnTotal*BnTotal*invRho*invRho));

  const CFreal ca = std::fabs(BnTotal)/sqrbar;
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
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
