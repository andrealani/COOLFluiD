#include "MHD2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

MHD2DVarSet::MHD2DVarSet(Common::SafePtr<BaseTerm> term) :
  ConvectiveVarSet(term),
  _model(term.d_castTo<MHDTerm>()),
  _BDipole(),
  _tanakaFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DVarSet::~MHD2DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DVarSet::computeFlux (const RealVector& data,
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
  const CFreal P = p + 0.5*B2;
  const CFreal E = p/(getModel()->getGamma() - 1.)
                   + 0.5 * (rho*(u*u + v*v) + B2);

  _fluxArray[0] = Vn*rho;
  _fluxArray[1] = Vn*rho*u - Bn*Bx + P*nx;
  _fluxArray[2] = Vn*rho*v - Bn*By + P*ny;
  _fluxArray[3] = - Bn*Bz;
  _fluxArray[4] = (v*Bx - By*u)*ny;
  _fluxArray[5] = (u*By - Bx*v)*nx;
  _fluxArray[6] = u*Bz*nx + v*Bz*ny;
  _fluxArray[7] = Vn*(E + P) - Bn*VdotB;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DVarSet::computeStateFlux (const RealVector& data)
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
  const CFreal P = p + 0.5*B2;
  const CFreal E = p/(getModel()->getGamma() - 1.)
                   + 0.5 * (rho*(u*u + v*v) + B2);
  const CFreal rhoU = rho*u;
  const CFreal rhoV = rho*v;
  const CFreal rhoUV = rhoU*v;
  const CFreal rhoUU = rhoU*u;
  const CFreal rhoVV = rhoV*v;

  _physFlux(0,XX) = rhoU;
  _physFlux(0,YY) = rhoV;

  _physFlux(1,XX) = rhoUU + Bx*Bx + P;
  _physFlux(1,YY) = rhoUV + BxBy;

  _physFlux(2,XX) = rhoUV + BxBy;
  _physFlux(2,YY) = rhoVV + By*By + P;

  _physFlux(3,XX) = -Bz*Bx;
  _physFlux(3,YY) = -Bz*By;

  _physFlux(4,XX) = 0;
  _physFlux(4,YY) = vBx - uBy;

  _physFlux(5,XX) = uBy - vBx;
  _physFlux(5,YY) = 0;

  _physFlux(6,XX) = u*Bz;
  _physFlux(6,YY) = v*Bz;

  _physFlux(7,XX) = u*(E + P) - Bx*VdotB;
  _physFlux(7,YY) = v*(E + P) - By*VdotB;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD2DVarSet::getMaxEigenValue(const RealVector& data,
				     const RealVector& normal)
{
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal B = data[MHDTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] +
    data[MHDTerm::VY]*normal[YY];
  const CFreal Bn = data[MHDTerm::BX]*normal[XX] +
    data[MHDTerm::BY]*normal[YY];
  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
                                  - 4.0*getModel()->getGamma()*p*Bn*Bn*
                                  invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  return Vn + cf;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD2DVarSet::getMaxAbsEigenValue(const RealVector& data,
					const RealVector& normal)
{
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal B = data[MHDTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] +
    data[MHDTerm::VY]*normal[YY];
  const CFreal Bn = data[MHDTerm::BX]*normal[XX] +
    data[MHDTerm::BY]*normal[YY];
  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
                                  - 4.0*getModel()->getGamma()*p*Bn*Bn*
                                  invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  return fabs(Vn) + cf;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DVarSet::computeMagneticDipole(CFreal xCoord, CFreal yCoord)
{
  // this formula is valid for a dipole at the origin
  const CFreal mX = getMX();
  const CFreal mY = getMY();
  const CFreal r = max(sqrt(xCoord*xCoord + yCoord*yCoord),(CFreal)1e-16); // this can fal in single precision
  const CFreal rCubeInv = 1.0/(r*r*r);
  const CFreal rSqInv = 1.0/(r*r);
  const CFreal mDotr = mX*xCoord + mY*yCoord;
  const CFreal BxDipole = rCubeInv*(3.0*xCoord*rSqInv*mDotr - mX);
  const CFreal ByDipole = rCubeInv*(3.0*yCoord*rSqInv*mDotr - mY);
  _BDipole[0] = BxDipole;
  _BDipole[1] = ByDipole;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DVarSet::computeEigenValues(const RealVector& data,
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

  result[0] = Vn - cf;
  result[1] = Vn - ca;
  result[2] = Vn - cs;
  result[3] = Vn;
  result[4] = Vn;
  result[5] = Vn + cs;
  result[6] = Vn + ca;
  result[7] = Vn + cf;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
