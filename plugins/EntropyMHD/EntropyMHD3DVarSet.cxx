#include "EntropyMHD3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DVarSet::EntropyMHD3DVarSet(Common::SafePtr<BaseTerm> term)
  : ConvectiveVarSet(term),
    _model(term.d_castTo<EntropyMHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DVarSet::~EntropyMHD3DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DVarSet::setup()
{
  ConvectiveVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DVarSet::computeFlux(const RealVector& pdata,
                                         const RealVector& normals)
{
  // Entropy-variable MHD + GLM flux projected onto normal direction.
  // Eqs 0-6,8: standard Dedner et al. (2002).
  // Eq 7: entropy advection sigma*Vn where sigma = p*rho^(1-gamma).
  // From d(sigma)/dt + div(sigma*v) = 0 (pure conservation law, no source term).
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];

  const CFreal rho   = pdata[EntropyMHDTerm::RHO];
  const CFreal u     = pdata[EntropyMHDTerm::VX];
  const CFreal v     = pdata[EntropyMHDTerm::VY];
  const CFreal w     = pdata[EntropyMHDTerm::VZ];
  const CFreal Bx    = pdata[EntropyMHDTerm::BX];
  const CFreal By    = pdata[EntropyMHDTerm::BY];
  const CFreal Bz    = pdata[EntropyMHDTerm::BZ];
  const CFreal p     = pdata[EntropyMHDTerm::P];
  const CFreal phi   = pdata[EntropyMHDTerm::PHI];
  const CFreal sigma = pdata[EntropyMHDTerm::T];  // sigma stored in T slot

  const CFreal Vn    = u*nx + v*ny + w*nz;
  const CFreal Bn    = Bx*nx + By*ny + Bz*nz;
  const CFreal sqB   = Bx*Bx + By*By + Bz*Bz;
  const CFreal Ptot  = p + 0.5*sqB;

  const CFreal refSpeedSq = _model->getRefSpeed() * _model->getRefSpeed();

  // Mass
  _fluxArray[0] = Vn*rho;

  // Momentum: rho*u_i*Vn - B_i*Bn + Ptot*n_i (uses recovered p in Ptot)
  _fluxArray[1] = Vn*rho*u - Bn*Bx + Ptot*nx;
  _fluxArray[2] = Vn*rho*v - Bn*By + Ptot*ny;
  _fluxArray[3] = Vn*rho*w - Bn*Bz + Ptot*nz;

  // Induction + GLM: (vB^T - Bv^T + psi*I) . n
  _fluxArray[4] = (Bx*v - By*u)*ny + (Bx*w - Bz*u)*nz + phi*nx;
  _fluxArray[5] = (By*u - Bx*v)*nx + (By*w - Bz*v)*nz + phi*ny;
  _fluxArray[6] = (Bz*u - Bx*w)*nx + (Bz*v - By*w)*ny + phi*nz;

  // Entropy advection: sigma * Vn (conservative, no source term)
  _fluxArray[7] = sigma * Vn;

  // GLM cleaning: ch^2 * Bn
  _fluxArray[8] = refSpeedSq * Bn;
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DVarSet::computeStateFlux(const RealVector& data)
{
  // Entropy-variable MHD + GLM 3D flux tensor.
  // Eq 7: sigma*v_i (entropy advection), eqs 0-6,8: standard MHD.
  const CFreal rho   = data[EntropyMHDTerm::RHO];
  const CFreal u     = data[EntropyMHDTerm::VX];
  const CFreal v     = data[EntropyMHDTerm::VY];
  const CFreal w     = data[EntropyMHDTerm::VZ];
  const CFreal Bx    = data[EntropyMHDTerm::BX];
  const CFreal By    = data[EntropyMHDTerm::BY];
  const CFreal Bz    = data[EntropyMHDTerm::BZ];
  const CFreal p     = data[EntropyMHDTerm::P];
  const CFreal phi   = data[EntropyMHDTerm::PHI];
  const CFreal sigma = data[EntropyMHDTerm::T];  // sigma stored in T slot

  const CFreal sqB   = Bx*Bx + By*By + Bz*Bz;
  const CFreal Ptot  = p + 0.5*sqB;
  const CFreal refSpeedSq = _model->getRefSpeed() * _model->getRefSpeed();

  // X-flux
  _physFlux(0, XX) = rho*u;
  _physFlux(1, XX) = rho*u*u - Bx*Bx + Ptot;
  _physFlux(2, XX) = rho*u*v - Bx*By;
  _physFlux(3, XX) = rho*u*w - Bx*Bz;
  _physFlux(4, XX) = phi;
  _physFlux(5, XX) = u*By - Bx*v;
  _physFlux(6, XX) = u*Bz - Bx*w;
  _physFlux(7, XX) = sigma * u;
  _physFlux(8, XX) = refSpeedSq * Bx;

  // Y-flux
  _physFlux(0, YY) = rho*v;
  _physFlux(1, YY) = rho*v*u - By*Bx;
  _physFlux(2, YY) = rho*v*v - By*By + Ptot;
  _physFlux(3, YY) = rho*v*w - By*Bz;
  _physFlux(4, YY) = v*Bx - By*u;
  _physFlux(5, YY) = phi;
  _physFlux(6, YY) = v*Bz - By*w;
  _physFlux(7, YY) = sigma * v;
  _physFlux(8, YY) = refSpeedSq * By;

  // Z-flux
  _physFlux(0, ZZ) = rho*w;
  _physFlux(1, ZZ) = rho*w*u - Bz*Bx;
  _physFlux(2, ZZ) = rho*w*v - Bz*By;
  _physFlux(3, ZZ) = rho*w*w - Bz*Bz + Ptot;
  _physFlux(4, ZZ) = w*Bx - Bz*u;
  _physFlux(5, ZZ) = w*By - Bz*v;
  _physFlux(6, ZZ) = phi;
  _physFlux(7, ZZ) = sigma * w;
  _physFlux(8, ZZ) = refSpeedSq * Bz;
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DVarSet::computeEigenValues(const RealVector& pdata,
                                                const RealVector& normal,
                                                RealVector& eValues)
{
  const CFreal rho   = pdata[EntropyMHDTerm::RHO];
  const CFreal Bx    = pdata[EntropyMHDTerm::BX];
  const CFreal By    = pdata[EntropyMHDTerm::BY];
  const CFreal Bz    = pdata[EntropyMHDTerm::BZ];
  const CFreal p     = pdata[EntropyMHDTerm::P];
  const CFreal gamma = pdata[EntropyMHDTerm::GAMMA];

  const CFreal Vn = pdata[EntropyMHDTerm::VX]*normal[XX]
                   + pdata[EntropyMHDTerm::VY]*normal[YY]
                   + pdata[EntropyMHDTerm::VZ]*normal[ZZ];

  const CFreal Bn  = Bx*normal[XX] + By*normal[YY] + Bz*normal[ZZ];
  const CFreal sqB = Bx*Bx + By*By + Bz*Bz;

  const CFreal astar2 = (gamma*p + sqB) / rho;
  const CFreal discrim = astar2*astar2 - 4.0*gamma*p*Bn*Bn / (rho*rho);
  const CFreal sqrtDiscrim = sqrt(std::max(discrim, 0.0));

  const CFreal cf = sqrt(0.5*(astar2 + sqrtDiscrim));
  const CFreal cs = sqrt(std::max(0.5*(astar2 - sqrtDiscrim), 0.0));
  const CFreal ca = std::abs(Bn) / sqrt(rho);

  const CFreal refSpeed = _model->getRefSpeed();

  eValues[0] = Vn - cf;
  eValues[1] = Vn - ca;
  eValues[2] = Vn - cs;
  eValues[3] = Vn;
  eValues[4] = Vn + cs;
  eValues[5] = Vn + ca;
  eValues[6] = Vn + cf;
  eValues[7] = refSpeed;
  eValues[8] = -refSpeed;
}

//////////////////////////////////////////////////////////////////////////////

CFreal EntropyMHD3DVarSet::getMaxEigenValue(const RealVector& pdata,
                                                const RealVector& normal)
{
  const CFreal rho   = pdata[EntropyMHDTerm::RHO];
  const CFreal Bx    = pdata[EntropyMHDTerm::BX];
  const CFreal By    = pdata[EntropyMHDTerm::BY];
  const CFreal Bz    = pdata[EntropyMHDTerm::BZ];
  const CFreal p     = pdata[EntropyMHDTerm::P];
  const CFreal gamma = pdata[EntropyMHDTerm::GAMMA];

  const CFreal Vn = pdata[EntropyMHDTerm::VX]*normal[XX]
                   + pdata[EntropyMHDTerm::VY]*normal[YY]
                   + pdata[EntropyMHDTerm::VZ]*normal[ZZ];

  const CFreal sqB = Bx*Bx + By*By + Bz*Bz;
  const CFreal Bn  = Bx*normal[XX] + By*normal[YY] + Bz*normal[ZZ];

  const CFreal astar2 = (gamma*p + sqB) / rho;
  const CFreal discrim = astar2*astar2 - 4.0*gamma*p*Bn*Bn / (rho*rho);
  const CFreal cf = sqrt(0.5*(astar2 + sqrt(std::max(discrim, 0.0))));

  const CFreal refSpeed = _model->getRefSpeed();

  return std::max(Vn + cf, refSpeed);
}

//////////////////////////////////////////////////////////////////////////////

CFreal EntropyMHD3DVarSet::getMaxAbsEigenValue(const RealVector& pdata,
                                                   const RealVector& normal)
{
  const CFreal rho   = pdata[EntropyMHDTerm::RHO];
  const CFreal Bx    = pdata[EntropyMHDTerm::BX];
  const CFreal By    = pdata[EntropyMHDTerm::BY];
  const CFreal Bz    = pdata[EntropyMHDTerm::BZ];
  const CFreal p     = pdata[EntropyMHDTerm::P];
  const CFreal gamma = pdata[EntropyMHDTerm::GAMMA];

  const CFreal Vn = pdata[EntropyMHDTerm::VX]*normal[XX]
                   + pdata[EntropyMHDTerm::VY]*normal[YY]
                   + pdata[EntropyMHDTerm::VZ]*normal[ZZ];

  const CFreal sqB = Bx*Bx + By*By + Bz*Bz;
  const CFreal Bn  = Bx*normal[XX] + By*normal[YY] + Bz*normal[ZZ];

  const CFreal astar2 = (gamma*p + sqB) / rho;
  const CFreal discrim = astar2*astar2 - 4.0*gamma*p*Bn*Bn / (rho*rho);
  const CFreal cf = sqrt(0.5*(astar2 + sqrt(std::max(discrim, 0.0))));

  const CFreal refSpeed = _model->getRefSpeed();

  return std::max(std::abs(Vn) + cf, refSpeed);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
