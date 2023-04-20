#include "MHD/MHD.hh"
#include "MHD2DProjectionPrim.hh"
#include "Common/CFLog.hh"
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

Environment::ObjectProvider<MHD2DProjectionPrim, ConvectiveVarSet, MHDModule, 1> 
mhd2DProjectionPrimProvider("MHD2DProjectionPrim");

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionPrim::MHD2DProjectionPrim(Common::SafePtr<BaseTerm> term) :
  MHD2DProjectionVarSet(term),
  _rm(),
  _rmInv(),
  _rightEv(),
  _leftEv()
{
  vector<std::string> names(9);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "p";
  names[8] = "phi";
  setVarNames(names);

  turnOnSourceTerm();
}

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionPrim::~MHD2DProjectionPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::setup()
{
  MHD2DProjectionVarSet::setup();

  _rm.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _rmInv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());

}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD2DProjectionPrim::getExtraVarNames() const
{
  vector<std::string> names(6);
  names[0] = "BxDipole";
  names[1] = "ByDipole";
  names[2] = "BxTotal";
  names[3] = "ByTotal";
  names[4] = "BTotal";
  names[5] = "rhoETotal";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::computeJacobians()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::computeEigenValuesVectors(RealMatrix& rightEv,
                                      RealMatrix& leftEv,
                                      RealVector& eValues,
                                      const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal rbar  = linearData[MHDProjectionTerm::RHO];
  const CFreal ubar  = linearData[MHDProjectionTerm::VX];
  const CFreal vbar  = linearData[MHDProjectionTerm::VY];
  const CFreal Bxbar = linearData[MHDProjectionTerm::BX];
  const CFreal Bybar = linearData[MHDProjectionTerm::BY];
  const CFreal Bzbar = linearData[MHDProjectionTerm::BZ];
  const CFreal pbar  = linearData[MHDProjectionTerm::P];
  const CFreal phibar = linearData[MHDProjectionTerm::PHI];
  const CFreal Vn = ubar*nx + vbar*ny;
  const CFreal Bn = Bxbar*nx + Bybar*ny;
  const CFreal Bp = Bybar*nx - Bxbar*ny;

  const CFreal gamma = getModel()->getGamma();
  const CFreal refSpeed = getModel()->getRefSpeed();

  CFLogDebugMax( "rbar = " << rbar << "\n");
  CFLogDebugMax( "ubar = " << ubar << "\n");
  CFLogDebugMax( "vbar = " << vbar << "\n");
  CFLogDebugMax( "Bxbar = " << Bxbar << "\n");
  CFLogDebugMax( "Bybar = " << Bybar << "\n");
  CFLogDebugMax( "Bzbar = " << Bzbar << "\n");
  CFLogDebugMax( "pbar = " << pbar << "\n");
  CFLogDebugMax( "phibar = " << phibar << "\n");
  CFLogDebugMax( "Vn = " << Vn << "\n");
  CFLogDebugMax( "Bn = " << Bn << "\n");
  CFLogDebugMax( "Bp = " << Bp << "\n");

  //  Characteristic speeds in the system

  const CFreal B2 = Bxbar*Bxbar + Bybar*Bybar + Bzbar*Bzbar;
  const CFreal astar2 = (gamma*pbar + B2)/rbar;

  CFLogDebugMax( "astar2 = " << astar2 << "\n");

  cf_assert(rbar > 0.);

  const CFreal sqrbar = sqrt(rbar);
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
          - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));
  CFreal cs2 = 0.5*(astar2 - sqrt(astar2*astar2
          - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));

  if (std::abs(cs2) < MathTools::MathConsts::CFrealEps()) cs2 = 0.;
  if (std::abs(cf2) < MathTools::MathConsts::CFrealEps()) cf2 = 0.;

  CFLogDebugMax( "cs2 = " << cs2 << "\n");
  CFLogDebugMax( "cf2 = " << cf2 << "\n");

  const CFreal cf = sqrt(cf2);
  const CFreal cs = sqrt(cs2);
  // const CFreal ca = Bn/sqrbar;
  const CFreal ca = std::abs(Bn)/sqrbar;
  const CFreal a  = linearData[MHDProjectionTerm::A];
  CFreal a2 = a*a;

  //fix against machine precision
  if (cs2 > a2) {
    a2 = cs2;
  }
  if (cf2 < a2) {
    a2 = cf2;
  }

  // Scaling factors for the normalized eigensystem

  CFreal alphaf2 = (a2 - cs2)/(cf2 - cs2);
  CFreal alphas2 = (cf2 - a2)/(cf2 - cs2);
  if (std::abs(alphaf2) < MathTools::MathConsts::CFrealEps()) alphaf2 = 0.;
  if (std::abs(alphas2) < MathTools::MathConsts::CFrealEps()) alphas2 = 0.;

  CFLogDebugMax( "alphas2 = " << alphas2 << "\n");
  CFLogDebugMax( "alphaf2 = " << alphaf2 << "\n");

  const CFreal alphaf = sqrt(alphaf2);
  const CFreal alphas = sqrt(alphas2);

  CFreal bethap = 0.;
  CFreal bethaz = 0.;

  if ((Bp*Bp + Bzbar*Bzbar) < MathTools::MathConsts::CFrealEps()) {
    bethap = 1.;
    bethaz = 0.;
  }
  else {
    bethap = Bp/sqrt(Bp*Bp + Bzbar*Bzbar);
    bethaz = Bzbar/sqrt(Bp*Bp + Bzbar*Bzbar);
  }

  const CFreal sgnbn = MathFunctions::sign(Bn);

  // sgnbn = 1. to make MHD work with B = 0

  const CFreal ovsq2 = 1./sqrt(2.);
  const CFreal ovsq2x = ovsq2*nx;
  const CFreal ovsq2y = ovsq2*ny;

  CFLogDebugMax( "sqrbar = " << sqrbar << "\n");
  CFLogDebugMax( "a = " << a << "\n");
  CFLogDebugMax( "bethap = " << bethap << "\n");
  CFLogDebugMax( "bethaz = " << bethaz << "\n");

  eValues[0] = Vn - cf;
  eValues[1] = Vn - ca;
  eValues[2] = Vn - cs;
  eValues[3] = -refSpeed;
  eValues[4] = Vn;
  eValues[5] = refSpeed;
  eValues[6] = Vn + cs;
  eValues[7] = Vn + ca;
  eValues[8] = Vn + cf;

  //       Left and right eigenmatrices
  // resize vector in the copy constructor of RealMatrix

  //  CFLogDebugMax( "R*L = " << "\n" << RealMatrix(rightEv*leftEv) << "\n");
  // CFLogDebugMax( "L*R = " << "\n" << RealMatrix(leftEv*rightEv) << "\n");

  const CFreal px = -ny;
  const CFreal py = nx;

  const CFreal asbppx = alphas*bethap*px;
  const CFreal afbppx = alphaf*bethap*px;
  const CFreal asbppy = alphas*bethap*py;
  const CFreal afbppy = alphaf*bethap*py;
  const CFreal asbz = alphas*cs*bethaz*sgnbn;
  const CFreal afbz = alphaf*cf*bethaz*sgnbn;
  const CFreal k1px = asbppx*cs*sgnbn;
  const CFreal k2px = afbppx*cf*sgnbn;
  const CFreal k1py = asbppy*cs*sgnbn;
  const CFreal k2py = afbppy*cf*sgnbn;
  const CFreal ov2a2 = 1./2./a2;
  const CFreal ascsnx = alphas*cs*nx;
  const CFreal afcfnx = alphaf*cf*nx;
  const CFreal ascsny = alphas*cs*ny;
  const CFreal afcfny = alphaf*cf*ny;
  const CFreal bpo = bethap*ovsq2;
  const CFreal bzo = bethaz*ovsq2;
  const CFreal bzopx = bzo*px;
  const CFreal bzopy = bzo*py;
  const CFreal sqra = sqrbar*a;
  const CFreal ov2sqra = 1./2./sqra;
  const CFreal ovsqrbar = 1./sqrbar;
  const CFreal bzsqra = sqra*bethaz;
  const CFreal ovbz = ov2sqra*bethaz;

  // matrix of right eigenvectors
  rightEv(0,0) = rbar*alphaf;
  rightEv(0,2) = rbar*alphas;
  rightEv(0,4) = 1.;
  rightEv(0,6) = rbar*alphas;
  rightEv(0,8) = rbar*alphaf;

  rightEv(1,0) = -afcfnx + k1px;
  rightEv(1,1) = -bzopx;
  rightEv(1,2) = -ascsnx - k2px;
  rightEv(1,6) = ascsnx + k2px;
  rightEv(1,7) = -bzopx;
  rightEv(1,8) = afcfnx - k1px;

  rightEv(2,0) = -afcfny + k1py;
  rightEv(2,1) = -bzopy;
  rightEv(2,2) = -ascsny - k2py;
  rightEv(2,6) = ascsny + k2py;
  rightEv(2,7) = -bzopy;
  rightEv(2,8) = afcfny - k1py;

  rightEv(3,0) = asbz;
  rightEv(3,1) = bpo;
  rightEv(3,2) = -afbz;
  rightEv(3,6) = afbz;
  rightEv(3,7) = bpo;
  rightEv(3,8) = -asbz;

  rightEv(4,0) = asbppx*sqra;
  rightEv(4,1) = -sqrbar*bzopx;
  rightEv(4,2) = -afbppx*sqra;
  rightEv(4,3) = ovsq2x;
  rightEv(4,5) = ovsq2x;
  rightEv(4,6) = -afbppx*sqra;
  rightEv(4,7) = sqrbar*bzopx;
  rightEv(4,8) = asbppx*sqra;

  rightEv(5,0) = asbppy*sqra;
  rightEv(5,1) = -sqrbar*bzopy;
  rightEv(5,2) = -afbppy*sqra;
  rightEv(5,3) = ovsq2y;
  rightEv(5,5) = ovsq2y;
  rightEv(5,6) = -afbppy*sqra;
  rightEv(5,7) = sqrbar*bzopy;
  rightEv(5,8) = asbppy*sqra;

  rightEv(6,0) = alphas*bzsqra;
  rightEv(6,1) = sqrbar*bpo;
  rightEv(6,2) = -alphaf*bzsqra;
  rightEv(6,6) = -alphaf*bzsqra;
  rightEv(6,7) = -sqrbar*bpo;
  rightEv(6,8) = alphas*bzsqra;

  rightEv(7,0) = alphaf*gamma*pbar;
  rightEv(7,2) = alphas*gamma*pbar;
  rightEv(7,6) = alphas*gamma*pbar;
  rightEv(7,8) = alphaf*gamma*pbar;

  rightEv(8,3) = -refSpeed*ovsq2;
  rightEv(8,5) = refSpeed*ovsq2;

  // matrix of left eigenvectors
  leftEv(0,1) = ov2a2*(-afcfnx + k1px);
  leftEv(0,2) = ov2a2*(-afcfny + k1py);
  leftEv(0,3) = ov2a2*asbz;
  leftEv(0,4) = asbppx*ov2sqra;
  leftEv(0,5) = asbppy*ov2sqra;
  leftEv(0,6) = alphas*ovbz;
  leftEv(0,7) = ov2a2*alphaf/rbar;

  leftEv(1,1) = -bzopx;
  leftEv(1,2) = -bzopy;
  leftEv(1,3) = bpo;
  leftEv(1,4) = -bzopx*ovsqrbar;
  leftEv(1,5) = -bzopy*ovsqrbar;
  leftEv(1,6) = bpo*ovsqrbar;

  leftEv(2,1) = ov2a2*(-ascsnx - k2px);
  leftEv(2,2) = ov2a2*(-ascsny - k2py);
  leftEv(2,3) = -ov2a2*afbz;
  leftEv(2,4) = -afbppx*ov2sqra;
  leftEv(2,5) = -afbppy*ov2sqra;
  leftEv(2,6) = -alphaf*ovbz;
  leftEv(2,7) = ov2a2*alphas/rbar;

  leftEv(3,4) = ovsq2x;
  leftEv(3,5) = ovsq2y;
  leftEv(3,8) = -1./refSpeed*ovsq2;

  leftEv(4,0) = 1.;
  leftEv(4,7) = -1./a2;

  leftEv(5,4) = ovsq2x;
  leftEv(5,5) = ovsq2y;
  leftEv(5,8) = 1./refSpeed*ovsq2;

  leftEv(6,1) = ov2a2*(ascsnx + k2px);
  leftEv(6,2) = ov2a2*(ascsny + k2py);
  leftEv(6,3) = ov2a2*afbz;
  leftEv(6,4) = -afbppx*ov2sqra;
  leftEv(6,5) = -afbppy*ov2sqra;
  leftEv(6,6) = -alphaf*ovbz;
  leftEv(6,7) = ov2a2*alphas/rbar;

  leftEv(7,1) = -bzopx;
  leftEv(7,2) = -bzopy;
  leftEv(7,3) = bpo;
  leftEv(7,4) = bzopx*ovsqrbar;
  leftEv(7,5) = bzopy*ovsqrbar;
  leftEv(7,6) = -bpo*ovsqrbar;

  leftEv(8,1) = ov2a2*(afcfnx -k1px);
  leftEv(8,2) = ov2a2*(afcfny -k1py);
  leftEv(8,3) = -ov2a2*asbz;
  leftEv(8,4) = asbppx*ov2sqra;
  leftEv(8,5) = asbppy*ov2sqra;
  leftEv(8,6) = alphas*ovbz;
  leftEv(8,7) = ov2a2*alphaf/rbar;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::setRotationMatrices(const RealVector& normal)
{
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal px = -ny;
  const CFreal py = nx;

  _rm(0,0) = 1.;
  _rm(1,1) = nx;
  _rm(1,2) = ny;
  _rm(2,1) = px;
  _rm(2,2) = py;
  _rm(3,3) = 1.;
  _rm(4,4) = nx;
  _rm(4,5) = ny;
  _rm(5,4) = px;
  _rm(5,5) = py;
  _rm(6,6) = 1.;
  _rm(7,7) = 1.;
  _rm(8,8) = 1.;


  _rmInv(0,0) = 1.;
  _rmInv(1,1) = nx;
  _rmInv(1,2) = px;
  _rmInv(2,1) = ny;
  _rmInv(2,2) = py;
  _rmInv(3,3) = 1.;
  _rmInv(4,4) = nx;
  _rmInv(4,5) = px;
  _rmInv(5,4) = ny;
  _rmInv(5,5) = py;
  _rmInv(6,6) = 1.;
  _rmInv(7,7) = 1.;
  _rmInv(8,8) = 1.;
}

//////////////////////////////////////////////////////////////////////////////

CFuint MHD2DProjectionPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::splitJacobian(RealMatrix& jacobPlus,
         RealMatrix& jacobMin,
         RealVector& eValues,
         const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal rbar  = linearData[MHDProjectionTerm::RHO];
  const CFreal ubar  = linearData[MHDProjectionTerm::VX];
  const CFreal vbar  = linearData[MHDProjectionTerm::VY];
  const CFreal Bxbar = linearData[MHDProjectionTerm::BX];
  const CFreal Bybar = linearData[MHDProjectionTerm::BY];
  const CFreal Bzbar = linearData[MHDProjectionTerm::BZ];
  const CFreal pbar  = linearData[MHDProjectionTerm::P];
  const CFreal phibar = linearData[MHDProjectionTerm::PHI];
  const CFreal Vn = ubar*nx + vbar*ny;
  const CFreal Bn = Bxbar*nx + Bybar*ny;
  const CFreal Bp = Bybar*nx - Bxbar*ny;
  const CFreal gamma = getModel()->getGamma();
  const CFreal refSpeed = getModel()->getRefSpeed();

  CFLogDebugMax( "rbar = " << rbar << "\n");
  CFLogDebugMax( "ubar = " << ubar << "\n");
  CFLogDebugMax( "vbar = " << vbar << "\n");
  CFLogDebugMax( "Bxbar = " << Bxbar << "\n");
  CFLogDebugMax( "Bybar = " << Bybar << "\n");
  CFLogDebugMax( "Bzbar = " << Bzbar << "\n");
  CFLogDebugMax( "pbar = " << pbar << "\n");
  CFLogDebugMax( "phibar = " << phibar << "\n");
  CFLogDebugMax( "Vn = " << Vn << "\n");
  CFLogDebugMax( "Bn = " << Bn << "\n");
  CFLogDebugMax( "Bp = " << Bp << "\n");

  //  Characteristic speeds in the system

  const CFreal B2 = Bxbar*Bxbar + Bybar*Bybar + Bzbar*Bzbar;
  const CFreal astar2 = (gamma*pbar + B2)/rbar;

  CFLogDebugMax( "astar2 = " << astar2 << "\n");

  cf_assert(rbar > 0.);

  const CFreal sqrbar = sqrt(rbar);
  CFreal cf2 = 0.5*(astar2 + sqrt(astar2*astar2
          - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));
  CFreal cs2 = 0.5*(astar2 - sqrt(astar2*astar2
          - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));

  if (std::abs(cs2) < MathTools::MathConsts::CFrealEps()) cs2 = 0.;
  if (std::abs(cf2) < MathTools::MathConsts::CFrealEps()) cf2 = 0.;

  CFLogDebugMax( "cs2 = " << cs2 << "\n");
  CFLogDebugMax( "cf2 = " << cs2 << "\n");

  const CFreal cf = sqrt(cf2);
  const CFreal cs = sqrt(cs2);

  // there was a bug in calculating the Alfven wave speed of sound
  // that should be tested

  // const CFreal ca = Bn/sqrbar;
  const CFreal ca = std::abs(Bn)/sqrbar;

  const CFreal a  = linearData[MHDProjectionTerm::A];
  CFreal a2 = a*a;

  //fix against machine precision
  if (cs2 > a2) {
    a2 = cs2;
  }
  if (cf2 < a2) {
    a2 = cf2;
  }

  // Scaling factors for the normalized eigensystem

  CFreal alphaf2 = (a2 - cs2)/(cf2 - cs2);
  CFreal alphas2 = (cf2 - a2)/(cf2 - cs2);
  if (std::abs(alphaf2) < MathTools::MathConsts::CFrealEps()) alphaf2 = 0.;
  if (std::abs(alphas2) < MathTools::MathConsts::CFrealEps()) alphas2 = 0.;

  CFLogDebugMax( "alphas2 = " << alphas2 << "\n");
  CFLogDebugMax( "alphaf2 = " << alphaf2 << "\n");

  const CFreal alphaf = sqrt(alphaf2);
  const CFreal alphas = sqrt(alphas2);

  CFreal bethap = 0.;
  CFreal bethaz = 0.;

  if ((Bp*Bp + Bzbar*Bzbar) < MathTools::MathConsts::CFrealEps()) {
    bethap = 1.;
    bethaz = 0.;
  }
  else {
    bethap = Bp/sqrt(Bp*Bp + Bzbar*Bzbar);
    bethaz = Bzbar/sqrt(Bp*Bp + Bzbar*Bzbar);
  }

  const CFreal sgnbn = MathFunctions::sign(Bn);

  // sgnbn = 1. to make MHD work with B = 0

  const CFreal ovsq2 = 1./sqrt(2.);
  const CFreal ovsq2x = ovsq2*nx;
  const CFreal ovsq2y = ovsq2*ny;

  CFLogDebugMax( "sqrbar = " << sqrbar << "\n");
  CFLogDebugMax( "a = " << a << "\n");
  CFLogDebugMax( "bethap = " << bethap << "\n");
  CFLogDebugMax( "bethaz = " << bethaz << "\n");

  eValues[0] = Vn - cf;
  eValues[1] = Vn - ca;
  eValues[2] = Vn - cs;
  eValues[3] = -refSpeed;
  eValues[4] = Vn;
  eValues[5] = refSpeed;
  eValues[6] = Vn + cs;
  eValues[7] = Vn + ca;
  eValues[8] = Vn + cf;

  //       Left and right eigenmatrices
  // resize vector in the copy constructor of RealMatrix

  //  CFLogDebugMax( "R*L = " << "\n" << RealMatrix(_rightEv*_leftEv) << "\n");
  // CFLogDebugMax( "L*R = " << "\n" << RealMatrix(_leftEv*_rightEv) << "\n");

  const CFreal px = -ny;
  const CFreal py = nx;

  const CFreal asbppx = alphas*bethap*px;
  const CFreal afbppx = alphaf*bethap*px;
  const CFreal asbppy = alphas*bethap*py;
  const CFreal afbppy = alphaf*bethap*py;
  const CFreal asbz = alphas*cs*bethaz*sgnbn;
  const CFreal afbz = alphaf*cf*bethaz*sgnbn;
  const CFreal k1px = asbppx*cs*sgnbn;
  const CFreal k2px = afbppx*cf*sgnbn;
  const CFreal k1py = asbppy*cs*sgnbn;
  const CFreal k2py = afbppy*cf*sgnbn;
  const CFreal ov2a2 = 1./2./a2;
  const CFreal ascsnx = alphas*cs*nx;
  const CFreal afcfnx = alphaf*cf*nx;
  const CFreal ascsny = alphas*cs*ny;
  const CFreal afcfny = alphaf*cf*ny;
  const CFreal bpo = bethap*ovsq2;
  const CFreal bzo = bethaz*ovsq2;
  const CFreal bzopx = bzo*px;
  const CFreal bzopy = bzo*py;
  const CFreal sqra = sqrbar*a;
  const CFreal ov2sqra = 1./2./sqra;
  const CFreal ovsqrbar = 1./sqrbar;
  const CFreal bzsqra = sqra*bethaz;
  const CFreal ovbz = ov2sqra*bethaz;

  // matrix of right eigenvectors
  _rightEv(0,0) = rbar*alphaf;
  _rightEv(0,2) = rbar*alphas;
  _rightEv(0,4) = 1.;
  _rightEv(0,6) = rbar*alphas;
  _rightEv(0,8) = rbar*alphaf;

  _rightEv(1,0) = -afcfnx + k1px;
  _rightEv(1,1) = -bzopx;
  _rightEv(1,2) = -ascsnx - k2px;
  _rightEv(1,6) = ascsnx + k2px;
  _rightEv(1,7) = -bzopx;
  _rightEv(1,8) = afcfnx - k1px;

  _rightEv(2,0) = -afcfny + k1py;
  _rightEv(2,1) = -bzopy;
  _rightEv(2,2) = -ascsny - k2py;
  _rightEv(2,6) = ascsny + k2py;
  _rightEv(2,7) = -bzopy;
  _rightEv(2,8) = afcfny - k1py;

  _rightEv(3,0) = asbz;
  _rightEv(3,1) = bpo;
  _rightEv(3,2) = -afbz;
  _rightEv(3,6) = afbz;
  _rightEv(3,7) = bpo;
  _rightEv(3,8) = -asbz;

  _rightEv(4,0) = asbppx*sqra;
  _rightEv(4,1) = -sqrbar*bzopx;
  _rightEv(4,2) = -afbppx*sqra;
  _rightEv(4,3) = ovsq2x;
  _rightEv(4,5) = ovsq2x;
  _rightEv(4,6) = -afbppx*sqra;
  _rightEv(4,7) = sqrbar*bzopx;
  _rightEv(4,8) = asbppx*sqra;

  _rightEv(5,0) = asbppy*sqra;
  _rightEv(5,1) = -sqrbar*bzopy;
  _rightEv(5,2) = -afbppy*sqra;
  _rightEv(5,3) = ovsq2y;
  _rightEv(5,5) = ovsq2y;
  _rightEv(5,6) = -afbppy*sqra;
  _rightEv(5,7) = sqrbar*bzopy;
  _rightEv(5,8) = asbppy*sqra;

  _rightEv(6,0) = alphas*bzsqra;
  _rightEv(6,1) = sqrbar*bpo;
  _rightEv(6,2) = -alphaf*bzsqra;
  _rightEv(6,6) = -alphaf*bzsqra;
  _rightEv(6,7) = -sqrbar*bpo;
  _rightEv(6,8) = alphas*bzsqra;

  _rightEv(7,0) = alphaf*gamma*pbar;
  _rightEv(7,2) = alphas*gamma*pbar;
  _rightEv(7,6) = alphas*gamma*pbar;
  _rightEv(7,8) = alphaf*gamma*pbar;

  _rightEv(8,3) = -refSpeed*ovsq2;
  _rightEv(8,5) = refSpeed*ovsq2;

  // matrix of left eigenvectors
  _leftEv(0,1) = ov2a2*(-afcfnx + k1px);
  _leftEv(0,2) = ov2a2*(-afcfny + k1py);
  _leftEv(0,3) = ov2a2*asbz;
  _leftEv(0,4) = asbppx*ov2sqra;
  _leftEv(0,5) = asbppy*ov2sqra;
  _leftEv(0,6) = alphas*ovbz;
  _leftEv(0,7) = ov2a2*alphaf/rbar;

  _leftEv(1,1) = -bzopx;
  _leftEv(1,2) = -bzopy;
  _leftEv(1,3) = bpo;
  _leftEv(1,4) = -bzopx*ovsqrbar;
  _leftEv(1,5) = -bzopy*ovsqrbar;
  _leftEv(1,6) = bpo*ovsqrbar;

  _leftEv(2,1) = ov2a2*(-ascsnx - k2px);
  _leftEv(2,2) = ov2a2*(-ascsny - k2py);
  _leftEv(2,3) = -ov2a2*afbz;
  _leftEv(2,4) = -afbppx*ov2sqra;
  _leftEv(2,5) = -afbppy*ov2sqra;
  _leftEv(2,6) = -alphaf*ovbz;
  _leftEv(2,7) = ov2a2*alphas/rbar;

  _leftEv(3,4) = ovsq2x;
  _leftEv(3,5) = ovsq2y;
  _leftEv(3,8) = -1./refSpeed*ovsq2;

  _leftEv(4,0) = 1.;
  _leftEv(4,7) = -1./a2;

  _leftEv(5,4) = ovsq2x;
  _leftEv(5,5) = ovsq2y;
  _leftEv(5,8) = 1./refSpeed*ovsq2;

  _leftEv(6,1) = ov2a2*(ascsnx + k2px);
  _leftEv(6,2) = ov2a2*(ascsny + k2py);
  _leftEv(6,3) = ov2a2*afbz;
  _leftEv(6,4) = -afbppx*ov2sqra;
  _leftEv(6,5) = -afbppy*ov2sqra;
  _leftEv(6,6) = -alphaf*ovbz;
  _leftEv(6,7) = ov2a2*alphas/rbar;

  _leftEv(7,1) = -bzopx;
  _leftEv(7,2) = -bzopy;
  _leftEv(7,3) = bpo;
  _leftEv(7,4) = bzopx*ovsqrbar;
  _leftEv(7,5) = bzopy*ovsqrbar;
  _leftEv(7,6) = -bpo*ovsqrbar;

  _leftEv(8,1) = ov2a2*(afcfnx -k1px);
  _leftEv(8,2) = ov2a2*(afcfny -k1py);
  _leftEv(8,3) = -ov2a2*asbz;
  _leftEv(8,4) = asbppx*ov2sqra;
  _leftEv(8,5) = asbppy*ov2sqra;
  _leftEv(8,6) = alphas*ovbz;
  _leftEv(8,7) = ov2a2*alphaf/rbar;

  // compute the eigen values + and -
  for (CFuint iEq = 0; iEq < 9; ++iEq) {
    _eValuesP[iEq] = max((CFreal)0.,eValues[iEq]);
    _eValuesM[iEq] = min((CFreal)0.,eValues[iEq]);
  }

  // compute jacobian + and -
  jacobPlus = _rightEv*(_eValuesP*_leftEv);
  jacobMin  = _rightEv*(_eValuesM*_leftEv);
 }

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::computePhysicalData(const State& state,
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
   data[MHDProjectionTerm::A] = sqrt(getModel()->getGamma()*
                                     data[MHDProjectionTerm::P]/rho);
  
   data[MHDProjectionTerm::GAMMA] = getModel()->getGamma();
 
   const RealVector& node = state.getCoordinates();
   data[MHDTerm::XP] =  node[XX];
   data[MHDTerm::YP] =  node[YY];
 }

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::setDimensionalValuesPlusExtraValues(const State& state,
                    RealVector& result,
                    RealVector& extra)
{
  RealVector BDipole(PhysicalModelStack::getActive()->getDim());

  const CFreal gammaMinus1 = getModel()->getGamma() - 1.;

  BDipole = getMagneticDipole(state.getCoordinates()[XX], state.getCoordinates()[YY]);
  const CFreal B1dotB0 = state[4]*BDipole[0] + state[5]*BDipole[1];
  const CFreal sqB0 = BDipole[0]*BDipole[0] + BDipole[1]*BDipole[1];
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal V2 = u*u + v*v;
  const CFreal sqB1 = state[4]*state[4] + state[5]*state[5];

  result[0] = state[0];
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
  result[4] = state[4];
  result[5] = state[5];
  result[6] = state[6];
  result[7] = state[7];
  result[8] = state[8];

  extra.resize(7);

  if ((getModel()->getMX() != 0.0) || (getModel()->getMY() != 0.0)) {
    extra[0] = BDipole[0];
    extra[1] = BDipole[1];
    extra[2] = state[4] + BDipole[0];
    extra[3] = state[5] + BDipole[1];
    extra[4] = sqrt(extra[2]*extra[2] + extra[3]*extra[3]);
    extra[5] = (state[7]/gammaMinus1) + 0.5*(rho*V2+sqB1) + B1dotB0 + 0.5*sqB0;
  }
  else {
    extra[0] = 0.0;
    extra[1] = 0.0;
    extra[2] = 0.0;
    extra[3] = 0.0;
    extra[4] = 0.0;
    extra[5] = 0.0;
  }

/*  std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName = nsp + "_divBNodal";

  DataHandle<CFreal> divBNodal = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

  const CFuint stateID = state.getLocalID();

  extra[6] = divBNodal[stateID];*/
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionPrim::computeStateFromPhysicalData(const RealVector& data,
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

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

