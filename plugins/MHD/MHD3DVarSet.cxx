#include "MHD3DVarSet.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

MHD3DVarSet::MHD3DVarSet(Common::SafePtr<Framework::BaseTerm> term) :
  ConvectiveVarSet(term),
  _model(term.d_castTo<MHDTerm>()),
  _BDipole(),
  _tanakaFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DVarSet::~MHD3DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DVarSet::computeFlux (const RealVector& data, const RealVector& normals)
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
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = u*nx  + v*ny + w*nz;
  const CFreal Bn = Bx*nx + By*ny + Bz*nz;
  const CFreal VdotB = u*Bx  + v*By  + w*Bz;
  const CFreal P = p + 0.5*(Bx*Bx + By*By + Bz*Bz);
  const CFreal E = p/(getModel()->getGamma() - 1.0) + 0.5*
    (rho*(u*u + v*v + w*w) + Bx*Bx + By*By + Bz*Bz);

  _fluxArray[0] = Vn*rho;
  _fluxArray[1] = Vn*rho*u - Bn*Bx + P*nx;
  _fluxArray[2] = Vn*rho*v - Bn*By + P*ny;
  _fluxArray[3] = Vn*rho*w - Bn*Bz + P*nz;
  _fluxArray[4] = (v*Bx - By*u)*ny + (w*Bx - Bz*u)*nz;
  _fluxArray[5] = (u*By - Bx*v)*nx + (w*By - Bz*v)*nz;
  _fluxArray[6] = (u*Bz - Bx*w)*nx + (v*Bz - By*w)*ny;
  _fluxArray[7] = Vn*(E + P) - Bn*VdotB;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DVarSet::computeStateFlux (const RealVector& data)
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
  const CFreal p = data[MHDTerm::P];
  const CFreal VdotB = u*Bx  + v*By  + w*Bz;
  const CFreal P = p + 0.5*(Bx*Bx + By*By + Bz*Bz);
  const CFreal E = p/(getModel()->getGamma() - 1.0) + 0.5*
    (rho*(u*u + v*v + w*w) + Bx*Bx + By*By + Bz*Bz);
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

  _physFlux(4,XX) = 0;
  _physFlux(4,YY) = vBx - uBy;
  _physFlux(4,ZZ) = wBx - uBz;

  _physFlux(5,XX) = uBy - vBx;
  _physFlux(5,YY) = 0;
  _physFlux(5,ZZ) = wBy - vBz;

  _physFlux(6,XX) = uBz - wBx;
  _physFlux(6,YY) = vBz - wBy;
  _physFlux(6,ZZ) = 0;

  _physFlux(7,XX) = u*(E + P) - Bx*VdotB;
  _physFlux(7,YY) = v*(E + P) - By*VdotB;
  _physFlux(7,ZZ) = w*(E + P) - Bz*VdotB;
}

//////////////////////////////////////////////////////////////////////

void MHD3DVarSet::setTransformationMatrices(const RealVector& coords,
			   RealVector& coordsSpherical,
                           RealMatrix& carSphTransMat,
                           RealMatrix& sphCarTransMat)
{
        const CFreal x = coords[0];
        const CFreal y = coords[1];
        const CFreal z = coords[2];
        const CFreal r = sqrt(x*x+y*y+z*z);
        const CFreal R = sqrt(x*x+y*y);

        coordsSpherical[0] = r;
        coordsSpherical[1] = atan2(R,z);
        coordsSpherical[2] = atan2(y,x);

        carSphTransMat(0,0) = x/r;
        carSphTransMat(1,0) = (x*z)/(R*r);
        carSphTransMat(2,0) = -y/R;

        carSphTransMat(0,1) = y/r;
        carSphTransMat(1,1) = (y*z)/(R*r);
        carSphTransMat(2,1) = x/R;

        carSphTransMat(0,2) = z/r;
        carSphTransMat(1,2) = -R/r;
        carSphTransMat(2,2) = 0.0;

        sphCarTransMat(0,0) = x/r;
        sphCarTransMat(1,0) = y/r;
        sphCarTransMat(2,0) = z/r;

        sphCarTransMat(0,1) = (x*z)/(R*r);
        sphCarTransMat(1,1) = (y*z)/(R*r);
        sphCarTransMat(2,1) = -R/r;

        sphCarTransMat(0,2) = -y/R;
        sphCarTransMat(1,2) = x/R;
        sphCarTransMat(2,2) = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHD3DVarSet::computeTanakaFluxPowell99Formulation(const RealVector& data,
                                                              const RealVector& normals)
{
  // According to the explanations in Powell, K.G. et. al., JCP, Vol.154, pp.284-309, 1999

  computeMagneticDipole(data[MHDTerm::XP],data[MHDTerm::YP],data[MHDTerm::ZP]);
  const CFreal Bx0 = _BDipole[0];
  const CFreal By0 = _BDipole[1];
  const CFreal Bz0 = _BDipole[2];

  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
  const CFreal rho = data[MHDTerm::RHO];
  const CFreal u = data[MHDTerm::VX];
  const CFreal v = data[MHDTerm::VY];
  const CFreal w = data[MHDTerm::VZ];
  const CFreal Bx1 = data[MHDTerm::BX];
  const CFreal By1 = data[MHDTerm::BY];
  const CFreal Bz1 = data[MHDTerm::BZ];

  // n stands for normal to the face of the control volume

  const CFreal Vn = u*nx + v*ny + w*nz;
  const CFreal Bn0 = Bx0*nx + By0*ny + Bz0*nz;

  const CFreal Bn1 = Bx1*nx + By1*ny + Bz1*nz;
  const CFreal sqB1 = Bx1*Bx1 + By1*By1 + Bz1*Bz1;
  const CFreal VdotB1 = u*Bx1  + v*By1 + w*Bz1;
  const CFreal B1dotB0 = Bx1*Bx0 + By1*By0 + Bz1*Bz0;
  const CFreal p = data[MHDTerm::P];

  const CFreal P = p + 0.5*sqB1 + B1dotB0;
  const CFreal E1 = p/(getModel()->getGamma() - 1.)
          + 0.5 * (rho*(u*u + v*v + w*w) + sqB1);

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _tanakaFlux[0] = Vn*rho;
  _tanakaFlux[1] = Vn*rho*u - Bn1*Bx1 - Bn0*Bx1 - Bn1*Bx0 + P*nx;
  _tanakaFlux[2] = Vn*rho*v - Bn1*By1 - Bn0*By1 - Bn1*By0 + P*ny;
  _tanakaFlux[3] = Vn*rho*w - Bn1*Bz1 - Bn0*Bz1 - Bn1*Bz0 + P*nz;
  _tanakaFlux[4] = (v*Bx1 - By1*u)*ny + (w*Bx1 - Bz1*u)*nz + (v*Bx0 - By0*u)*ny + (w*Bx0 - Bz0*u)*nz;
  _tanakaFlux[5] = (u*By1 - Bx1*v)*nx + (w*By1 - Bz1*v)*nz + (u*By0 - Bx0*v)*nx + (w*By0 - Bz0*v)*nz;
  _tanakaFlux[6] = (u*Bz1 - Bx1*w)*nx + (v*Bz1 - By1*w)*ny + (u*Bz0 - Bx0*w)*nx + (v*Bz0 - By0*w)*ny;
  _tanakaFlux[7] = Vn*(E1 + p + 0.5*sqB1) - Bn1*VdotB1 + Vn*B1dotB0 - Bn0*VdotB1;

  return _tanakaFlux;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DVarSet::getMaxEigenValue(const RealVector& data,
				     const RealVector& normal)
{
  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal B = data[MHDTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] +
    data[MHDTerm::VY]*normal[YY] + data[MHDTerm::VZ]*normal[ZZ];

  const CFreal Bn = data[MHDTerm::BX]*normal[XX] +
    data[MHDTerm::BY]*normal[YY] + data[MHDTerm::BZ]*normal[ZZ];

  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
				  - 4.0*getModel()->getGamma()*p*Bn*Bn*
				  invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  return Vn + cf;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DVarSet::getMaxAbsEigenValue(const RealVector& data,
					const RealVector& normal)
{  const CFreal invRho = 1./data[MHDTerm::RHO];
  const CFreal B = data[MHDTerm::B];
  const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] +
    data[MHDTerm::VY]*normal[YY] + data[MHDTerm::VZ]*normal[ZZ];

  const CFreal Bn = data[MHDTerm::BX]*normal[XX] +
    data[MHDTerm::BY]*normal[YY] + data[MHDTerm::BZ]*normal[ZZ];

  const CFreal astar2 = (getModel()->getGamma()*p + B2)*invRho;
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
          - 4.0*getModel()->getGamma()*p*Bn*Bn*
          invRho*invRho));

  const CFreal cf = sqrt(std::abs(cf2));
  return fabs(Vn) + cf;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DVarSet::computeMagneticDipole(CFreal xCoord, CFreal yCoord, CFreal zCoord)
{
  // this formula is valid for a dipole at the origin
  const CFreal mX = getMX();
  const CFreal mY = getMY();
  const CFreal mZ = getMZ();

  const CFreal r = sqrt(xCoord*xCoord + yCoord*yCoord + zCoord*zCoord);
  const CFreal rCubeInv = 1.0/(r*r*r);
  const CFreal rSqInv = 1.0/(r*r);
  const CFreal mDotr = mX*xCoord + mY*yCoord + mZ*zCoord;

  const CFreal BxDipole = rCubeInv*(3.0*xCoord*rSqInv*mDotr - mX);
  const CFreal ByDipole = rCubeInv*(3.0*yCoord*rSqInv*mDotr - mY);
  const CFreal BzDipole = rCubeInv*(3.0*zCoord*rSqInv*mDotr - mZ);

  _BDipole[0] = BxDipole;
  _BDipole[1] = ByDipole;
  _BDipole[2] = BzDipole;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DVarSet::computeEigenValues(const RealVector& data,
				     const RealVector& normal,
				     RealVector& result)
{
  // eigenvalues based on total magnetic field, Btotal
  computeMagneticDipole(data[MHDTerm::XP],data[MHDTerm::YP],data[MHDTerm::ZP]);
  const CFreal Bx0 = _BDipole[0];
  const CFreal By0 = _BDipole[1];
  const CFreal Bz0 = _BDipole[2];
  const CFreal Bx1 = data[MHDTerm::BX];
  const CFreal By1 = data[MHDTerm::BY];
  const CFreal Bz1 = data[MHDTerm::BZ];
  const CFreal BxTotal = Bx1 + Bx0;
  const CFreal ByTotal = By1 + By0;
  const CFreal BzTotal = Bz1 + Bz0;

  const CFreal invRho = 1./data[MHDTerm::RHO];
  //const CFreal B = data[MHDTerm::B];
  const CFreal sqBTotal = BxTotal*BxTotal + ByTotal*ByTotal + BzTotal*BzTotal;
  //const CFreal B2 = B*B;
  const CFreal p = data[MHDTerm::P];
  const CFreal Vn = data[MHDTerm::VX]*normal[XX] + data[MHDTerm::VY]*normal[YY] +
    data[MHDTerm::VZ]*normal[ZZ];
  //const CFreal Bn = data[MHDTerm::BX]*normal[XX] + data[MHDTerm::BY]*normal[YY] +
  //data[MHDTerm::BZ]*normal[ZZ];
  const CFreal BnTotal = BxTotal*normal[XX] + ByTotal*normal[YY] + BzTotal*normal[ZZ];
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
