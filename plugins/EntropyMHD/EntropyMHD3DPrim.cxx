#include "EntropyMHD/EntropyMHD.hh"
#include "EntropyMHD3DPrim.hh"
#include "Common/CFLog.hh"
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

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EntropyMHD3DPrim, ConvectiveVarSet, EntropyMHDModule, 1>
entropyMHD3DPrimProvider("EntropyMHD3DPrim");

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DPrim::EntropyMHD3DPrim(Common::SafePtr<BaseTerm> term)
  : EntropyMHD3DVarSet(term)
{
  vector<std::string> names(9);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "sigma";
  names[8] = "phi";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DPrim::~EntropyMHD3DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DPrim::setup()
{
  EntropyMHD3DVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DPrim::computePhysicalData(const State& state,
                                     RealVector& data)
{
  const CFreal rho   = state[0];
  const CFreal u     = state[1];
  const CFreal v     = state[2];
  const CFreal w     = state[3];
  const CFreal Bx    = state[4];
  const CFreal By    = state[5];
  const CFreal Bz    = state[6];
  const CFreal sigma = state[7];  // entropy variable sigma = p*rho^(1-gamma)
  const CFreal phi   = state[8];

  const CFreal gamma = getModel()->getGamma();

  // Recover thermodynamic pressure: p = sigma * rho^(gamma-1)
  const CFreal p = sigma * std::pow(rho, gamma - 1.0);

  data[EntropyMHDTerm::RHO]   = rho;
  data[EntropyMHDTerm::VX]    = u;
  data[EntropyMHDTerm::VY]    = v;
  data[EntropyMHDTerm::VZ]    = w;
  data[EntropyMHDTerm::BX]    = Bx;
  data[EntropyMHDTerm::BY]    = By;
  data[EntropyMHDTerm::BZ]    = Bz;
  data[EntropyMHDTerm::P]     = p;
  data[EntropyMHDTerm::PHI]   = phi;
  data[EntropyMHDTerm::GAMMA] = gamma;
  data[EntropyMHDTerm::T]     = sigma;  // sigma slot (named T for historical reasons)

  // Coordinates (read by EntropyMHDSourceTerm for gravity)
  const RealVector& coord = state.getCoordinates();
  data[EntropyMHDTerm::XP] = coord[XX];
  data[EntropyMHDTerm::YP] = coord[YY];
  data[EntropyMHDTerm::ZP] = coord[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DPrim::computeStateFromPhysicalData(const RealVector& data,
                                              State& state)
{
  state[0] = data[EntropyMHDTerm::RHO];
  state[1] = data[EntropyMHDTerm::VX];
  state[2] = data[EntropyMHDTerm::VY];
  state[3] = data[EntropyMHDTerm::VZ];
  state[4] = data[EntropyMHDTerm::BX];
  state[5] = data[EntropyMHDTerm::BY];
  state[6] = data[EntropyMHDTerm::BZ];
  // state[7] = sigma = p * rho^(1-gamma)
  state[7] = data[EntropyMHDTerm::P] * std::pow(data[EntropyMHDTerm::RHO], 1.0 - data[EntropyMHDTerm::GAMMA]);
  state[8] = data[EntropyMHDTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  // Public entry point: linearize about the model's stored pdata.
  computeEigenValuesVectorsFromPData(getModel()->getPhysicalData(),
                                      rightEv, leftEv, eValues, normal);
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DPrim::computeEigenValuesVectorsFromPData(const RealVector& linearData,
                                                    RealMatrix& rightEv,
                                                    RealMatrix& leftEv,
                                                    RealVector& eValues,
                                                    const RealVector& normal)
{
  // Strategy:
  //   1. Build the 9-wave Powell-Dedner GLM eigensystem in MHD primitive
  //      basis (..., p, ...) and MHD wave ordering, using the same closed-form
  //      formulas as MHD3DProjectionPrim::computeEigenValuesVectors. The
  //      pressure pbar comes from EntropyMHD pdata[P] = sigma * rho^(gamma-1) (set
  //      by computePhysicalData), so the wave speeds are physically the same
  //      as in MHD.
  //   2. Apply the sigma-vs-p basis transformation T (see class header):
  //         R_EntropyMHD = T^{-1} * R_MHD  -> only row 7 changes
  //         L_EntropyMHD = L_MHD * T       -> only columns 0 and 7 change
  //   3. Permute MHD wave ordering -> EntropyMHD wave ordering
  //         MHD:        [Vn-cf, Vn-ca, Vn-cs, -c_h, Vn,   +c_h, Vn+cs, Vn+ca, Vn+cf]
  //         EntropyMHD: [Vn-cf, Vn-ca, Vn-cs,  Vn, Vn+cs, Vn+ca, Vn+cf, +c_h, -c_h ]

  // ---- Step 0: tangent basis ------------------------------------------------

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  CFreal px = 0., py = 0., pz = 0.;
  CFreal zx = 0., zy = 0., zz = 0.;
  CFreal length_x_y = std::sqrt(nx*nx + ny*ny);
  CFreal over_length_x_y = 0.;

  if (length_x_y < MathTools::MathConsts::CFrealEps()) {
    // Normal nearly aligned with y-axis: rotate the (p, z) tangent frame.
    length_x_y = std::sqrt(nx*nx + nz*nz);
    over_length_x_y = 1.0/length_x_y;

    px =  over_length_x_y*nz;
    py =  0.0;
    pz = -over_length_x_y*nx;

    zx = -over_length_x_y*ny*nx;
    zy =  length_x_y;
    zz = -over_length_x_y*ny*nz;
  }
  else {
    over_length_x_y = 1.0/length_x_y;

    px =  over_length_x_y*ny;
    py = -over_length_x_y*nx;
    pz =  0.0;

    zx =  over_length_x_y*nx*nz;
    zy =  over_length_x_y*ny*nz;
    zz = -length_x_y;
  }

  const CFreal rbar  = linearData[EntropyMHDTerm::RHO];
  const CFreal ubar  = linearData[EntropyMHDTerm::VX];
  const CFreal vbar  = linearData[EntropyMHDTerm::VY];
  const CFreal wbar  = linearData[EntropyMHDTerm::VZ];
  const CFreal Bxbar = linearData[EntropyMHDTerm::BX];
  const CFreal Bybar = linearData[EntropyMHDTerm::BY];
  const CFreal Bzbar = linearData[EntropyMHDTerm::BZ];
  const CFreal pbar  = linearData[EntropyMHDTerm::P];   // p = sigma * rho^(gamma-1); set in computePhysicalData
  const CFreal sbar  = linearData[EntropyMHDTerm::T];   // sigma; direct entropy variable
  // const CFreal phibar = linearData[EntropyMHDTerm::PHI];  // unused: phi only enters via the c_h waves

  const CFreal Vn = ubar*nx + vbar*ny + wbar*nz;
  const CFreal Bn = Bxbar*nx + Bybar*ny + Bzbar*nz;
  const CFreal Bp = Bxbar*px + Bybar*py + Bzbar*pz;
  const CFreal Bz_ = Bxbar*zx + Bybar*zy + Bzbar*zz;  // renamed to avoid shadowing zx,zy,zz components

  const CFreal gamma = getModel()->getGamma();
  const CFreal refSpeed = getModel()->getRefSpeed();

  cf_assert(rbar > 0.);
  cf_assert(pbar > 0.);

  // ---- Step 1: characteristic speeds (cf, ca, cs, sound a) ------------------
  // (Identical math to MHD3DProjectionPrim; only the source of pbar differs.)
  const CFreal B2 = Bxbar*Bxbar + Bybar*Bybar + Bzbar*Bzbar;
  const CFreal astar2 = (gamma*pbar + B2)/rbar;

  const CFreal sqrbar = std::sqrt(rbar);
  // Bit-identical to EntropyMHD3DVarSet::computeEigenValues so that the eigenvalues
  // returned here match `computeEigenValues` exactly under FP precision
  // (matters at PIL collapse where Bn -> 0 and cs -> 0 with very small
  // sub-EPS magnitude; guarding cs2 with std::abs<EPS would zero it,
  // disagreeing with the reference by ~|cs| at the same level).
  const CFreal discrimEntropyMHD = astar2*astar2 - 4.0*gamma*pbar*Bn*Bn/(rbar*rbar);
  const CFreal sqrtDiscrim = std::sqrt(std::max(0.0, discrimEntropyMHD));
  CFreal cf2 = 0.5*(astar2 + sqrtDiscrim);
  CFreal cs2 = std::max(0.5*(astar2 - sqrtDiscrim), 0.0);

  const CFreal cf = std::sqrt(cf2);
  const CFreal cs = std::sqrt(cs2);
  const CFreal ca = std::abs(Bn)/sqrbar;

  // Sound speed from pressure (a^2 = gamma*p/rho).
  CFreal a2 = gamma*pbar/rbar;
  CFreal a  = std::sqrt(a2);

  // Machine-precision fixes (identical to MHD).
  if (cs2 > a2) a2 = cs2;
  if (cf2 < a2) a2 = cf2;

  // ---- Step 2: wave-mode scaling factors -----------------------------------
  CFreal alphaf2 = (cf2 - cs2 != 0.0) ? (a2 - cs2)/(cf2 - cs2) : 1.0;
  CFreal alphas2 = (cf2 - cs2 != 0.0) ? (cf2 - a2)/(cf2 - cs2) : 0.0;
  if (std::abs(alphaf2) < MathTools::MathConsts::CFrealEps()) alphaf2 = 0.;
  if (std::abs(alphas2) < MathTools::MathConsts::CFrealEps()) alphas2 = 0.;

  const CFreal alphaf = std::sqrt(alphaf2);
  const CFreal alphas = std::sqrt(alphas2);

  CFreal bethap = 1.;
  CFreal bethaz = 0.;
  if ((Bp*Bp + Bz_*Bz_) >= MathTools::MathConsts::CFrealEps()) {
    const CFreal invBpz = 1.0/std::sqrt(Bp*Bp + Bz_*Bz_);
    bethap = Bp*invBpz;
    bethaz = Bz_*invBpz;
  }

  const CFreal sgnbn = MathFunctions::sign(Bn);

  const CFreal ovsq2  = 1.0/std::sqrt(2.0);
  const CFreal ovsq2x = ovsq2*nx;
  const CFreal ovsq2y = ovsq2*ny;
  const CFreal ovsq2z = ovsq2*nz;

  // ---- Step 3: MHD-ordered eigenvalues -------------------------------------
  // (We will permute to EntropyMHD order at the very end.)
  RealVector eVals_MHD(9);
  eVals_MHD[0] = Vn - cf;
  eVals_MHD[1] = Vn - ca;
  eVals_MHD[2] = Vn - cs;
  eVals_MHD[3] = -refSpeed;
  eVals_MHD[4] = Vn;
  eVals_MHD[5] = refSpeed;
  eVals_MHD[6] = Vn + cs;
  eVals_MHD[7] = Vn + ca;
  eVals_MHD[8] = Vn + cf;

  // ---- Step 4: shorthand combinations (verbatim from MHD3DProjectionPrim) --
  const CFreal asbppx = alphas*bethap*px;
  const CFreal afbppx = alphaf*bethap*px;
  const CFreal asbppy = alphas*bethap*py;
  const CFreal afbppy = alphaf*bethap*py;
  const CFreal asbppz = alphas*bethap*pz;
  const CFreal afbppz = alphaf*bethap*pz;
  const CFreal asbzzx = alphas*bethaz*zx;
  const CFreal afbzzx = alphaf*bethaz*zx;
  const CFreal asbzzy = alphas*bethaz*zy;
  const CFreal afbzzy = alphaf*bethaz*zy;
  const CFreal asbzzz = alphas*bethaz*zz;
  const CFreal afbzzz = alphaf*bethaz*zz;
  const CFreal asbz   = alphas*cs*bethaz*sgnbn;  // unused below; kept for parity with MHD source
  (void)asbz;
  const CFreal k1px = asbppx*cs*sgnbn;
  const CFreal k2px = afbppx*cf*sgnbn;
  const CFreal k1py = asbppy*cs*sgnbn;
  const CFreal k2py = afbppy*cf*sgnbn;
  const CFreal k1pz = asbppz*cs*sgnbn;
  const CFreal k2pz = afbppz*cf*sgnbn;
  const CFreal k1zx = asbzzx*cs*sgnbn;
  const CFreal k2zx = afbzzx*cf*sgnbn;
  const CFreal k1zy = asbzzy*cs*sgnbn;
  const CFreal k2zy = afbzzy*cf*sgnbn;
  const CFreal k1zz = asbzzz*cs*sgnbn;
  const CFreal k2zz = afbzzz*cf*sgnbn;

  const CFreal ov2a2 = (a2 > 0.0) ? 1.0/(2.0*a2) : 0.0;
  const CFreal ascsnx = alphas*cs*nx;
  const CFreal afcfnx = alphaf*cf*nx;
  const CFreal ascsny = alphas*cs*ny;
  const CFreal afcfny = alphaf*cf*ny;
  const CFreal ascsnz = alphas*cs*nz;
  const CFreal afcfnz = alphaf*cf*nz;
  const CFreal bpo = bethap*ovsq2;
  const CFreal bzo = bethaz*ovsq2;
  const CFreal bzopx = bzo*px;
  const CFreal bzopy = bzo*py;
  const CFreal bzopz = bzo*pz;
  const CFreal bpozx = bpo*zx;
  const CFreal bpozy = bpo*zy;
  const CFreal bpozz = bpo*zz;
  const CFreal sqra = sqrbar*a;
  const CFreal ov2sqra = (sqra > 0.0) ? 1.0/(2.0*sqra) : 0.0;
  const CFreal ovsqrbar = 1.0/sqrbar;

  // ---- Step 5: right eigenvectors in MHD ordering and p-basis --------------
  RealMatrix rEv_MHD(9, 9, 0.0);

  // col 0 (Vn - cf)
  rEv_MHD(0,0) = rbar*alphaf;
  rEv_MHD(1,0) = -afcfnx + k1px + k1zx;
  rEv_MHD(2,0) = -afcfny + k1py + k1zy;
  rEv_MHD(3,0) = -afcfnz + k1pz + k1zz;
  rEv_MHD(4,0) =  asbppx*sqra + asbzzx*sqra;
  rEv_MHD(5,0) =  asbppy*sqra + asbzzy*sqra;
  rEv_MHD(6,0) =  asbppz*sqra + asbzzz*sqra;
  rEv_MHD(7,0) =  alphaf*gamma*pbar;

  // col 1 (Vn - ca)
  rEv_MHD(1,1) = -bzopx + bethap*ovsq2*zx;
  rEv_MHD(2,1) = -bzopy + bethap*ovsq2*zy;
  rEv_MHD(3,1) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  rEv_MHD(4,1) = -sqrbar*bzopx + sqrbar*ovsq2*bethap*zx;
  rEv_MHD(5,1) = -sqrbar*bzopy + sqrbar*ovsq2*bethap*zy;
  rEv_MHD(6,1) = -sqrbar*bzopz + sqrbar*ovsq2*bethap*zz;

  // col 2 (Vn - cs)
  rEv_MHD(0,2) =  rbar*alphas;
  rEv_MHD(1,2) = -ascsnx - k2px - k2zx;
  rEv_MHD(2,2) = -ascsny - k2py - k2zy;
  rEv_MHD(3,2) = -ascsnz - k2pz - k2zz;
  rEv_MHD(4,2) = -afbppx*sqra - afbzzx*sqra;
  rEv_MHD(5,2) = -afbppy*sqra - afbzzy*sqra;
  rEv_MHD(6,2) = -afbppz*sqra - afbzzz*sqra;
  rEv_MHD(7,2) =  alphas*gamma*pbar;

  // col 3 (-c_h)
  rEv_MHD(4,3) = ovsq2x;
  rEv_MHD(5,3) = ovsq2y;
  rEv_MHD(6,3) = ovsq2z;
  rEv_MHD(8,3) = -refSpeed*ovsq2;

  // col 4 (Vn; entropy/contact)
  rEv_MHD(0,4) = 1.0;

  // col 5 (+c_h)
  rEv_MHD(4,5) = ovsq2x;
  rEv_MHD(5,5) = ovsq2y;
  rEv_MHD(6,5) = ovsq2z;
  rEv_MHD(8,5) =  refSpeed*ovsq2;

  // col 6 (Vn + cs)
  rEv_MHD(0,6) =  rbar*alphas;
  rEv_MHD(1,6) =  ascsnx + k2px + k2zx;
  rEv_MHD(2,6) =  ascsny + k2py + k2zy;
  rEv_MHD(3,6) =  ascsnz + k2pz + k2zz;
  rEv_MHD(4,6) = -afbppx*sqra - afbzzx*sqra;
  rEv_MHD(5,6) = -afbppy*sqra - afbzzy*sqra;
  rEv_MHD(6,6) = -afbppz*sqra - afbzzz*sqra;
  rEv_MHD(7,6) =  alphas*gamma*pbar;

  // col 7 (Vn + ca)
  rEv_MHD(1,7) = -bzopx + bethap*ovsq2*zx;
  rEv_MHD(2,7) = -bzopy + bethap*ovsq2*zy;
  rEv_MHD(3,7) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  rEv_MHD(4,7) =  sqrbar*bzopx - sqrbar*ovsq2*bethap*zx;
  rEv_MHD(5,7) =  sqrbar*bzopy - sqrbar*ovsq2*bethap*zy;
  rEv_MHD(6,7) =  sqrbar*bzopz - sqrbar*ovsq2*bethap*zz;

  // col 8 (Vn + cf)
  rEv_MHD(0,8) =  rbar*alphaf;
  rEv_MHD(1,8) =  afcfnx - k1px - k1zx;
  rEv_MHD(2,8) =  afcfny - k1py - k1zy;
  rEv_MHD(3,8) =  alphaf*cf*nz - k1pz - k1zz;
  rEv_MHD(4,8) =  asbppx*sqra + asbzzx*sqra;
  rEv_MHD(5,8) =  asbppy*sqra + asbzzy*sqra;
  rEv_MHD(6,8) =  asbppz*sqra + asbzzz*sqra;
  rEv_MHD(7,8) =  alphaf*gamma*pbar;

  // ---- Step 6: left eigenvectors in MHD ordering and p-basis ---------------
  RealMatrix lEv_MHD(9, 9, 0.0);

  // row 0 (Vn - cf)
  lEv_MHD(0,1) = ov2a2*(-afcfnx + k1px + k1zx);
  lEv_MHD(0,2) = ov2a2*(-afcfny + k1py + k1zy);
  lEv_MHD(0,3) = ov2a2*(-afcfnz + k1pz + k1zz);
  lEv_MHD(0,4) = ov2sqra*(asbppx + asbzzx);
  lEv_MHD(0,5) = ov2sqra*(asbppy + asbzzy);
  lEv_MHD(0,6) = ov2sqra*(asbppz + asbzzz);
  lEv_MHD(0,7) = ov2a2*alphaf/rbar;

  // row 1 (Vn - ca)
  lEv_MHD(1,1) = -bzopx + bpozx;
  lEv_MHD(1,2) = -bzopy + bpozy;
  lEv_MHD(1,3) = -bzopz + bpozz;
  lEv_MHD(1,4) = ovsqrbar*(-bzopx + bpozx);
  lEv_MHD(1,5) = ovsqrbar*(-bzopy + bpozy);
  lEv_MHD(1,6) = ovsqrbar*(-bzopz + bpozz);

  // row 2 (Vn - cs)
  lEv_MHD(2,1) = ov2a2*(-ascsnx - k2px - k2zx);
  lEv_MHD(2,2) = ov2a2*(-ascsny - k2py - k2zy);
  lEv_MHD(2,3) = ov2a2*(-ascsnz - k2pz - k2zz);
  lEv_MHD(2,4) = ov2sqra*(-afbppx - afbzzx);
  lEv_MHD(2,5) = ov2sqra*(-afbppy - afbzzy);
  lEv_MHD(2,6) = ov2sqra*(-afbppz - afbzzz);
  lEv_MHD(2,7) = ov2a2*alphas/rbar;

  // row 3 (-c_h)
  lEv_MHD(3,4) =  ovsq2x;
  lEv_MHD(3,5) =  ovsq2y;
  lEv_MHD(3,6) =  ovsq2z;
  lEv_MHD(3,8) = -ovsq2/refSpeed;

  // row 4 (Vn; entropy/contact)
  lEv_MHD(4,0) = 1.0;
  lEv_MHD(4,7) = (a2 > 0.0) ? -1.0/a2 : 0.0;

  // row 5 (+c_h)
  lEv_MHD(5,4) = ovsq2x;
  lEv_MHD(5,5) = ovsq2y;
  lEv_MHD(5,6) = ovsq2z;
  lEv_MHD(5,8) = ovsq2/refSpeed;

  // row 6 (Vn + cs)
  lEv_MHD(6,1) = ov2a2*(ascsnx + k2px + k2zx);
  lEv_MHD(6,2) = ov2a2*(ascsny + k2py + k2zy);
  lEv_MHD(6,3) = ov2a2*(ascsnz + k2pz + k2zz);
  lEv_MHD(6,4) = ov2sqra*(-afbppx - afbzzx);
  lEv_MHD(6,5) = ov2sqra*(-afbppy - afbzzy);
  lEv_MHD(6,6) = ov2sqra*(-afbppz - afbzzz);
  lEv_MHD(6,7) = ov2a2*alphas/rbar;

  // row 7 (Vn + ca)
  lEv_MHD(7,1) = -bzopx + bpozx;
  lEv_MHD(7,2) = -bzopy + bpozy;
  lEv_MHD(7,3) = -bzopz + bpozz;
  lEv_MHD(7,4) = ovsqrbar*( bzopx - bpozx);
  lEv_MHD(7,5) = ovsqrbar*( bzopy - bpozy);
  lEv_MHD(7,6) = ovsqrbar*( bzopz - bpozz);

  // row 8 (Vn + cf)
  lEv_MHD(8,1) = ov2a2*( afcfnx - k1px - k1zx);
  lEv_MHD(8,2) = ov2a2*( afcfny - k1py - k1zy);
  lEv_MHD(8,3) = ov2a2*( afcfnz - k1pz - k1zz);
  lEv_MHD(8,4) = ov2sqra*(asbppx + asbzzx);
  lEv_MHD(8,5) = ov2sqra*(asbppy + asbzzy);
  lEv_MHD(8,6) = ov2sqra*(asbppz + asbzzz);
  lEv_MHD(8,7) = ov2a2*alphaf/rbar;

  // ---- Step 7: sigma-vs-p basis transformation ------------------------------
  // R_EntropyMHD = T^{-1} * R_MHD : modify row 7 only
  //   R_EntropyMHD[7,k] = - (gamma-1) * sigma / rho * R_MHD[0,k]
  //                       + rho^(1-gamma)            * R_MHD[7,k]
  const CFreal coefRrho   = -(gamma - 1.0) * sbar / rbar;
  const CFreal coefRsigma =  std::pow(rbar, 1.0 - gamma);
  for (CFuint k = 0; k < 9; ++k) {
    const CFreal R0k = rEv_MHD(0, k);
    const CFreal R7k = rEv_MHD(7, k);
    rEv_MHD(7, k) = coefRrho * R0k + coefRsigma * R7k;
  }

  // L_EntropyMHD = L_MHD * T : modify columns 0 (rho) and 7 (sigma) only
  //   L_EntropyMHD[k,0] = L_MHD[k,0] + (gamma-1) * p / rho * L_MHD[k,7]
  //   L_EntropyMHD[k,7] = rho^(gamma-1) * L_MHD[k,7]
  const CFreal coefLrho   = (gamma - 1.0) * pbar / rbar;
  const CFreal coefLsigma = std::pow(rbar, gamma - 1.0);
  for (CFuint k = 0; k < 9; ++k) {
    const CFreal L_k_p = lEv_MHD(k, 7);
    lEv_MHD(k, 0) += coefLrho * L_k_p;
    lEv_MHD(k, 7)  = coefLsigma * L_k_p;
  }

  // ---- Step 8: permute MHD ordering -> EntropyMHD ordering -----------------
  // entropyMHD_to_mhd[k_EntropyMHD] gives the MHD column/row index for EntropyMHD slot k_EntropyMHD.
  static const CFuint entropyMHD_to_mhd[9] = {0, 1, 2, 4, 6, 7, 8, 5, 3};

  for (CFuint k = 0; k < 9; ++k) {
    const CFuint km = entropyMHD_to_mhd[k];
    eValues[k] = eVals_MHD[km];
    for (CFuint i = 0; i < 9; ++i) {
      rightEv(i, k) = rEv_MHD(i, km);
    }
    for (CFuint j = 0; j < 9; ++j) {
      leftEv(k, j) = lEv_MHD(km, j);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
