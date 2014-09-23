#include <numeric>

#include "MathTools/MathFunctions.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "MHD/MHD.hh"
#include "MHD/MHD3DPrim.hh"
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

Environment::ObjectProvider<MHD3DPrim, ConvectiveVarSet, MHDModule, 1> mhd3DPrimProvider("MHD3DPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DPrim::MHD3DPrim(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DVarSet(term),
  _tempR(),
  _tempL(),
  _rm(),
  _rmInv(),
  _rightEv(),
  _leftEv()
{
  vector<std::string> names(8);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "p";
  setVarNames(names);

  turnOnSourceTerm();
}

//////////////////////////////////////////////////////////////////////////////

MHD3DPrim::~MHD3DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::setup()
{
  MHD3DVarSet::setup();

  _tempR.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _tempL.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _rm.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _rmInv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());

}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD3DPrim::getExtraVarNames() const
{
  vector<std::string> names(9);
  names[0] = "BxDipole";
  names[1] = "ByDipole";
  names[2] = "BzDipole";
  names[3] = "BxTotal";
  names[4] = "ByTotal";
  names[5] = "BzTotal";
  names[6] = "BTotal";
  names[7] = "rhoETotal";
  names[8] = "divB";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::computeJacobians()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
              RealMatrix& leftEv,
              RealVector& eValues,
              const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  CFreal px = 0.;
  CFreal py = 0.;
  CFreal pz = 0.;
  CFreal zx = 0.;
  CFreal zy = 0.;
  CFreal zz = 0.;
  CFreal length_x_y = sqrt(nx*nx + ny*ny);
  CFreal over_length_x_y = 0.;

  if (length_x_y < MathTools::MathConsts::CFrealEps()) {
    length_x_y = sqrt(nx*nx + nz*nz);
    over_length_x_y = 1./length_x_y;

    px =  over_length_x_y*nz;
    py = 0.;
    pz = -over_length_x_y*nx;

    zx =  -over_length_x_y*ny*nx;
    zy =  length_x_y;
    zz = -over_length_x_y*ny*nz;
  }
  else {
    over_length_x_y = 1./length_x_y;

    px =  over_length_x_y*ny;
    py = -over_length_x_y*nx;
    pz = 0.;

    zx =  over_length_x_y*nx*nz;
    zy =  over_length_x_y*ny*nz;
    zz = -length_x_y;
  }

  const CFreal gamma = getModel()->getGamma();
  const CFreal rbar  = linearData[MHDTerm::RHO];
  const CFreal ubar  = linearData[MHDTerm::VX];
  const CFreal vbar  = linearData[MHDTerm::VY];
  const CFreal wbar  = linearData[MHDTerm::VZ];
  const CFreal Bxbar = linearData[MHDTerm::BX];
  const CFreal Bybar = linearData[MHDTerm::BY];
  const CFreal Bzbar = linearData[MHDTerm::BZ];
  const CFreal pbar  = linearData[MHDTerm::P];

  const CFreal Vn = ubar*nx  + vbar*ny + wbar*nz;
  // const CFreal Vp = ubar*px  + vbar*py + wbar*pz;
  // const CFreal Vz = ubar*zx  + vbar*zy + wbar*zz;

  const CFreal Bn = Bxbar*nx + Bybar*ny + Bzbar*nz;;
  const CFreal Bp = Bxbar*px + Bybar*py + Bzbar*pz;
  const CFreal Bz = Bxbar*zx + Bybar*zy + Bzbar*zz;

  //  Characteristic speeds in the system

  const CFreal B2 = Bxbar*Bxbar + Bybar*Bybar + Bzbar*Bzbar;
  const CFreal astar2 = (gamma*pbar + B2)/rbar;
  const CFreal sqrbar = sqrt(rbar);
  CFreal cf2 = 0.5* (astar2 + sqrt(astar2*astar2
           - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));
  CFreal cs2 = 0.5* (astar2 - sqrt(astar2*astar2
           - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));

  if (std::abs(cs2) < MathTools::MathConsts::CFrealEps()) cs2 = 0.;
  if (std::abs(cf2) < MathTools::MathConsts::CFrealEps()) cf2 = 0.;

  CFLogDebugMax( "cs2 = " << cs2 << "\n");
  CFLogDebugMax( "cf2 = " << cs2 << "\n");

  const CFreal cf = sqrt(cf2);
  const CFreal cs = sqrt(cs2);

  // there was a bug in calculating the Alfven wave speed of sound that should be tested

  //const CFreal ca = Bn/sqrbar;
  const CFreal ca = std::abs(Bn)/sqrbar;

  CFreal a  = linearData[MHDTerm::A];
  CFreal a2 = a*a;

  //fix against machine precision
  if (cs2 > a2) {
    a2 = cs2;
  }
  if (cf2 < a2) {
    a2 = cf2;
  }
  // Scaling factors for the normaized eigensystem

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

  if ((Bp*Bp + Bz*Bz) < MathTools::MathConsts::CFrealEps()) {
    bethap = 1.;
    bethaz = 0.;
  }
  else {
    bethap = Bp/sqrt(Bp*Bp + Bz*Bz);
    bethaz = Bz/sqrt(Bp*Bp + Bz*Bz);
  }

  const CFreal sgnbn = MathFunctions::sign(Bn);
  const CFreal ovsq2 = 1./sqrt(2.);

  eValues[0] = Vn - cf;
  eValues[1] = Vn - ca;
  eValues[2] = Vn - cs;
  eValues[3] = Vn;
  eValues[4] = Vn;
  eValues[5] = Vn + cs;
  eValues[6] = Vn + ca;
  eValues[7] = Vn + cf;

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
  const CFreal asbz = alphas*cs*bethaz*sgnbn;
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

  const CFreal ov2a2 = 1./2./a2;
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
  const CFreal ov2sqra = 1./2./sqra;
  const CFreal ovsqrbar = 1./sqrbar;

  // matrix of right eigenvectors
  rightEv(0,0) = rbar*alphaf;
  rightEv(0,2) = rbar*alphas;
  rightEv(0,3) = 1.;
  rightEv(0,5) = rbar*alphas;
  rightEv(0,7) = rbar*alphaf;

  rightEv(1,0) = -afcfnx + k1px + k1zx;
  rightEv(1,1) = -bzopx + bethap*ovsq2*zx;
  rightEv(1,2) = -ascsnx - k2px - k2zx;
  rightEv(1,5) = ascsnx + k2px + k2zx;
  rightEv(1,6) = -bzopx + bethap*ovsq2*zx;
  rightEv(1,7) = afcfnx - k1px - k1zx;

  rightEv(2,0) = -afcfny + k1py + k1zy;
  rightEv(2,1) = -bzopy + bethap*ovsq2*zy;
  rightEv(2,2) = -ascsny - k2py - k2zy;
  rightEv(2,5) = ascsny + k2py + k2zy;
  rightEv(2,6) = -bzopy + bethap*ovsq2*zy;
  rightEv(2,7) = afcfny - k1py - k1zy;

  rightEv(3,0) = -afcfnz + k1pz + asbz*zz;
  rightEv(3,1) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  rightEv(3,2) = -alphas*cs*nz - k2pz - k2zz;
  rightEv(3,5) = alphas*cs*nz + k2pz + k2zz;
  rightEv(3,6) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  rightEv(3,7) = alphaf*cf*nz - k1pz - k1zz;

  rightEv(4,0) = asbppx*sqra + asbzzx*sqra;
  rightEv(4,1) = -sqrbar*bzopx + sqrbar*ovsq2*bethap*zx;
  rightEv(4,2) = -afbppx*sqra - afbzzx*sqra;
  rightEv(4,4) = nx;
  rightEv(4,5) = -afbppx*sqra - afbzzx*sqra;
  rightEv(4,6) = sqrbar*bzopx - sqrbar*ovsq2*bethap*zx;
  rightEv(4,7) = asbppx*sqra + asbzzx*sqra;

  rightEv(5,0) = asbppy*sqra + asbzzy*sqra;
  rightEv(5,1) = -sqrbar*bzopy + sqrbar*ovsq2*bethap*zy;
  rightEv(5,2) = -afbppy*sqra - afbzzy*sqra;
  rightEv(5,4) = ny;
  rightEv(5,5) = -afbppy*sqra - afbzzy*sqra;
  rightEv(5,6) = sqrbar*bzopy - sqrbar*ovsq2*bethap*zy;
  rightEv(5,7) = asbppy*sqra + asbzzy*sqra;

  rightEv(6,0) = asbppz*sqra + asbzzz*sqra;
  rightEv(6,1) = -sqrbar*bzopz + sqrbar*ovsq2*bethap*zz;
  rightEv(6,2) = -afbppz*sqra - afbzzz*sqra;
  rightEv(6,4) = nz;
  rightEv(6,5) = -afbppz*sqra - afbzzz*sqra;
  rightEv(6,6) = sqrbar*bzopz - sqrbar*ovsq2*bethap*zz;
  rightEv(6,7) = asbppz*sqra + asbzzz*sqra;

  rightEv(7,0) = alphaf*gamma*pbar;
  rightEv(7,2) = alphas*gamma*pbar;
  rightEv(7,5) = alphas*gamma*pbar;
  rightEv(7,7) = alphaf*gamma*pbar;

  // matrix of left eigenvectors
  leftEv(0,1) = ov2a2*(-afcfnx + k1px + k1zx);
  leftEv(0,2) = ov2a2*(-afcfny + k1py + k1zy);
  leftEv(0,3) = ov2a2*(-afcfnz + k1pz + k1zz);
  leftEv(0,4) = ov2sqra*(asbppx + asbzzx);
  leftEv(0,5) = ov2sqra*(asbppy + asbzzy);
  leftEv(0,6) = ov2sqra*(asbppz + asbzzz);
  leftEv(0,7) = ov2a2*alphaf/rbar;

  leftEv(1,1) = -bzopx + bpozx;
  leftEv(1,2) = -bzopy + bpozy;
  leftEv(1,3) = -bzopz + bpozz;
  leftEv(1,4) = ovsqrbar*(-bzopx + bpozx);
  leftEv(1,5) = ovsqrbar*(-bzopy + bpozy);
  leftEv(1,6) = ovsqrbar*(-bzopz + bpozz);

  leftEv(2,1) = ov2a2*(-ascsnx - k2px - k2zx);
  leftEv(2,2) = ov2a2*(-ascsny - k2py - k2zy);
  leftEv(2,3) = ov2a2*(-ascsnz - k2pz - k2zz);
  leftEv(2,4) = ov2sqra*(-afbppx - afbzzx);
  leftEv(2,5) = ov2sqra*(-afbppy - afbzzy);
  leftEv(2,6) = ov2sqra*(-afbppz - afbzzz);
  leftEv(2,7) = ov2a2*alphas/rbar;

  leftEv(3,0) = 1.;
  leftEv(3,7) = -1./a2;

  leftEv(4,4) = nx;
  leftEv(4,5) = ny;
  leftEv(4,6) = nz;

  leftEv(5,1) = ov2a2*(ascsnx + k2px + k2zx);
  leftEv(5,2) = ov2a2*(ascsny + k2py + k2zy);
  leftEv(5,3) = ov2a2*(ascsnz + k2pz + k2zz);
  leftEv(5,4) = ov2sqra*(-afbppx - afbzzx);
  leftEv(5,5) = ov2sqra*(-afbppy - afbzzy);
  leftEv(5,6) = ov2sqra*(-afbppz - afbzzz);
  leftEv(5,7) = ov2a2*alphas/rbar;

  leftEv(6,1) = -bzopx + bpozx;
  leftEv(6,2) = -bzopy + bpozy;
  leftEv(6,3) = -bzopz + bpozz;
  leftEv(6,4) = ovsqrbar*(bzopx - bpozx);
  leftEv(6,5) = ovsqrbar*(bzopy - bpozy);
  leftEv(6,6) = ovsqrbar*(bzopz - bpozz);

  leftEv(7,1) = ov2a2*(afcfnx - k1px - k1zx);
  leftEv(7,2) = ov2a2*(afcfny - k1py - k1zy);
  leftEv(7,3) = ov2a2*(afcfnz - k1pz - k1zz);
  leftEv(7,4) = ov2sqra*(asbppx + asbzzx);
  leftEv(7,5) = ov2sqra*(asbppy + asbzzy);
  leftEv(7,6) = ov2sqra*(asbppz + asbzzz);
  leftEv(7,7) = ov2a2*alphaf/rbar;

  //  CFLogDebugMax( "R*L = " << "\n" << RealMatrix(rightEv*leftEv) << "\n");
  // CFLogDebugMax( "L*R = " << "\n" << RealMatrix(leftEv*rightEv) << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::setRotationMatrices(const RealVector& normals)
{
  const CFreal nx = normals[0];
  const CFreal ny = normals[1];
  const CFreal nz = normals[2];

  CFreal px = 0.;
  CFreal py = 0.;
  CFreal pz = 0.;
  CFreal zx = 0.;
  CFreal zy = 0.;
  CFreal zz = 0.;
  CFreal length_x_y = sqrt(nx*nx + ny*ny);
  CFreal over_length_x_y = 0.;

  if (length_x_y < MathTools::MathConsts::CFrealEps()) {
    length_x_y = sqrt(nx*nx + nz*nz);
    over_length_x_y = 1./length_x_y;

    px =  over_length_x_y*nz;
    py = 0.;
    pz = -over_length_x_y*nx;

    zx =  -over_length_x_y*ny*nx;
    zy =  length_x_y;
    zz = -over_length_x_y*ny*nz;
  }
  else {
    over_length_x_y = 1./length_x_y;

    px =  over_length_x_y*ny;
    py = -over_length_x_y*nx;
    pz = 0.;

    zx =  over_length_x_y*nx*nz;
    zy =  over_length_x_y*ny*nz;
    zz = -length_x_y;
  }

  _rm(0,0) = 1.;

  _rm(1,1) = nx;
  _rm(1,2) = ny;
  _rm(1,3) = nz;

  _rm(2,1) = px;
  _rm(2,2) = py;
  _rm(2,3) = pz;

  _rm(3,1) = zx;
  _rm(3,2) = zy;
  _rm(3,3) = zz;

  _rm(4,4) = nx;
  _rm(4,5) = ny;
  _rm(4,6) = nz;

  _rm(5,4) = px;
  _rm(5,5) = py;
  _rm(5,6) = pz;

  _rm(6,4) = zx;
  _rm(6,5) = zy;
  _rm(6,6) = zz;

  _rm(7,7) = 1.;

  _rmInv(0,0) = 1.;

  _rmInv(1,1) = nx;
  _rmInv(1,2) = px;
  _rmInv(1,3) = zx;

  _rmInv(2,1) = ny;
  _rmInv(2,2) = py;
  _rmInv(2,3) = zy;

  _rmInv(3,1) = nz;
  _rmInv(3,2) = pz;
  _rmInv(3,3) = zz;

  _rmInv(4,4) = nx;
  _rmInv(4,5) = px;
  _rmInv(4,6) = zx;

  _rmInv(5,4) = ny;
  _rmInv(5,5) = py;
  _rmInv(5,6) = zy;

  _rmInv(6,4) = nz;
  _rmInv(6,5) = pz;
  _rmInv(6,6) = zz;

  _rmInv(7,7) = 1.;

  //  CFLogDebugMax( "RM*RMI = " << "\n" << RealMatrix(_rm*_rmInv) << "\n");
  // CFLogDebugMax( "RMI*RM = " << "\n" << RealMatrix(_rmInv*_rm) << "\n");
}

//////////////////////////////////////////////////////////////////////////////

CFuint MHD3DPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::splitJacobian(RealMatrix& jacobPlus,
         RealMatrix& jacobMin,
         RealVector& eValues,
         const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[0];
  const CFreal ny = normal[1];
  const CFreal nz = normal[2];
  CFreal px = 0.;
  CFreal py = 0.;
  CFreal pz = 0.;
  CFreal zx = 0.;
  CFreal zy = 0.;
  CFreal zz = 0.;
  CFreal length_x_y = sqrt(nx*nx + ny*ny);
  CFreal over_length_x_y = 0.;

  if (length_x_y < MathTools::MathConsts::CFrealEps()) {
    length_x_y = sqrt(nx*nx + nz*nz);
    over_length_x_y = 1./length_x_y;

    px =  over_length_x_y*nz;
    py = 0.;
    pz = -over_length_x_y*nx;

    zx =  -over_length_x_y*ny*nx;
    zy =  length_x_y;
    zz = -over_length_x_y*ny*nz;
  }
  else {
    over_length_x_y = 1./length_x_y;

    px =  over_length_x_y*ny;
    py = -over_length_x_y*nx;
    pz = 0.;

    zx =  over_length_x_y*nx*nz;
    zy =  over_length_x_y*ny*nz;
    zz = -length_x_y;
  }

  const CFreal gamma = getModel()->getGamma();
  const CFreal rbar  = linearData[MHDTerm::RHO];
  const CFreal ubar  = linearData[MHDTerm::VX];
  const CFreal vbar  = linearData[MHDTerm::VY];
  const CFreal wbar  = linearData[MHDTerm::VZ];
  const CFreal Bxbar = linearData[MHDTerm::BX];
  const CFreal Bybar = linearData[MHDTerm::BY];
  const CFreal Bzbar = linearData[MHDTerm::BZ];
  const CFreal pbar  = linearData[MHDTerm::P];
  const CFreal Vn = ubar*nx  + vbar*ny + wbar*nz;
  // const CFreal Vp = ubar*px  + vbar*py + wbar*pz;
  // const CFreal Vz = ubar*zx  + vbar*zy + wbar*zz;

  const CFreal Bn = Bxbar*nx + Bybar*ny + Bzbar*nz;;
  const CFreal Bp = Bxbar*px + Bybar*py + Bzbar*pz;
  const CFreal Bz = Bxbar*zx + Bybar*zy + Bzbar*zz;

  //  Characteristic speeds in the system

  const CFreal B2 = Bxbar*Bxbar + Bybar*Bybar + Bzbar*Bzbar;
  const CFreal astar2 = (gamma*pbar + B2)/rbar;
  const CFreal sqrbar = sqrt(rbar);
  CFreal cf2 = 0.5* (astar2 + sqrt(astar2*astar2
           - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));
  CFreal cs2 = 0.5* (astar2 - sqrt(astar2*astar2
           - 4.0*gamma*pbar*Bn*Bn/rbar/rbar));

  if (std::abs(cs2) < MathTools::MathConsts::CFrealEps()) cs2 = 0.;
  if (std::abs(cf2) < MathTools::MathConsts::CFrealEps()) cf2 = 0.;

  CFLogDebugMax( "cs2 = " << cs2 << "\n");
  CFLogDebugMax( "cf2 = " << cs2 << "\n");

  const CFreal cf = sqrt(cf2);
  const CFreal cs = sqrt(cs2);
  const CFreal ca = std::abs(Bn)/sqrbar;
  CFreal a  = linearData[MHDTerm::A];
  CFreal a2 = a*a;

  //fix against machine precision
  if (cs2 > a2) {
    a2 = cs2;
  }
  if (cf2 < a2) {
    a2 = cf2;
  }
  // Scaling factors for the normaized eigensystem

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

  if ((Bp*Bp + Bz*Bz) < MathTools::MathConsts::CFrealEps()) {
    bethap = 1.;
    bethaz = 0.;
  }
  else {
    bethap = Bp/sqrt(Bp*Bp + Bz*Bz);
    bethaz = Bz/sqrt(Bp*Bp + Bz*Bz);
  }

  const CFreal sgnbn = MathFunctions::sign(Bn);
  const CFreal ovsq2 = 1./sqrt(2.);

  eValues[0] = Vn - cf;
  eValues[1] = Vn - ca;
  eValues[2] = Vn - cs;
  eValues[3] = Vn;
  eValues[4] = Vn;
  eValues[5] = Vn + cs;
  eValues[6] = Vn + ca;
  eValues[7] = Vn + cf;

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
  const CFreal asbz = alphas*cs*bethaz*sgnbn;
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

  const CFreal ov2a2 = 1./2./a2;
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
  const CFreal ov2sqra = 1./2./sqra;
  const CFreal ovsqrbar = 1./sqrbar;

  // matrix of right eigenvectors
  _rightEv(0,0) = rbar*alphaf;
  _rightEv(0,2) = rbar*alphas;
  _rightEv(0,3) = 1.;
  _rightEv(0,5) = rbar*alphas;
  _rightEv(0,7) = rbar*alphaf;

  _rightEv(1,0) = -afcfnx + k1px + k1zx;
  _rightEv(1,1) = -bzopx + bethap*ovsq2*zx;
  _rightEv(1,2) = -ascsnx - k2px - k2zx;
  _rightEv(1,5) = ascsnx + k2px + k2zx;
  _rightEv(1,6) = -bzopx + bethap*ovsq2*zx;
  _rightEv(1,7) = afcfnx - k1px - k1zx;

  _rightEv(2,0) = -afcfny + k1py + k1zy;
  _rightEv(2,1) = -bzopy + bethap*ovsq2*zy;
  _rightEv(2,2) = -ascsny - k2py - k2zy;
  _rightEv(2,5) = ascsny + k2py + k2zy;
  _rightEv(2,6) = -bzopy + bethap*ovsq2*zy;
  _rightEv(2,7) = afcfny - k1py - k1zy;

  _rightEv(3,0) = -afcfnz + k1pz + asbz*zz;
  _rightEv(3,1) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  _rightEv(3,2) = -alphas*cs*nz - k2pz - k2zz;
  _rightEv(3,5) = alphas*cs*nz + k2pz + k2zz;
  _rightEv(3,6) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  _rightEv(3,7) = alphaf*cf*nz - k1pz - k1zz;

  _rightEv(4,0) = asbppx*sqra + asbzzx*sqra;
  _rightEv(4,1) = -sqrbar*bzopx + sqrbar*ovsq2*bethap*zx;
  _rightEv(4,2) = -afbppx*sqra - afbzzx*sqra;
  _rightEv(4,4) = nx;
  _rightEv(4,5) = -afbppx*sqra - afbzzx*sqra;
  _rightEv(4,6) = sqrbar*bzopx - sqrbar*ovsq2*bethap*zx;
  _rightEv(4,7) = asbppx*sqra + asbzzx*sqra;

  _rightEv(5,0) = asbppy*sqra + asbzzy*sqra;
  _rightEv(5,1) = -sqrbar*bzopy + sqrbar*ovsq2*bethap*zy;
  _rightEv(5,2) = -afbppy*sqra - afbzzy*sqra;
  _rightEv(5,4) = ny;
  _rightEv(5,5) = -afbppy*sqra - afbzzy*sqra;
  _rightEv(5,6) = sqrbar*bzopy - sqrbar*ovsq2*bethap*zy;
  _rightEv(5,7) = asbppy*sqra + asbzzy*sqra;

  _rightEv(6,0) = asbppz*sqra + asbzzz*sqra;
  _rightEv(6,1) = -sqrbar*bzopz + sqrbar*ovsq2*bethap*zz;
  _rightEv(6,2) = -afbppz*sqra - afbzzz*sqra;
  _rightEv(6,4) = nz;
  _rightEv(6,5) = -afbppz*sqra - afbzzz*sqra;
  _rightEv(6,6) = sqrbar*bzopz - sqrbar*ovsq2*bethap*zz;
  _rightEv(6,7) = asbppz*sqra + asbzzz*sqra;

  _rightEv(7,0) = alphaf*gamma*pbar;
  _rightEv(7,2) = alphas*gamma*pbar;
  _rightEv(7,5) = alphas*gamma*pbar;
  _rightEv(7,7) = alphaf*gamma*pbar;

  // matrix of left eigenvectors
  _leftEv(0,1) = ov2a2*(-afcfnx + k1px + k1zx);
  _leftEv(0,2) = ov2a2*(-afcfny + k1py + k1zy);
  _leftEv(0,3) = ov2a2*(-afcfnz + k1pz + k1zz);
  _leftEv(0,4) = ov2sqra*(asbppx + asbzzx);
  _leftEv(0,5) = ov2sqra*(asbppy + asbzzy);
  _leftEv(0,6) = ov2sqra*(asbppz + asbzzz);
  _leftEv(0,7) = ov2a2*alphaf/rbar;

  _leftEv(1,1) = -bzopx + bpozx;
  _leftEv(1,2) = -bzopy + bpozy;
  _leftEv(1,3) = -bzopz + bpozz;
  _leftEv(1,4) = ovsqrbar*(-bzopx + bpozx);
  _leftEv(1,5) = ovsqrbar*(-bzopy + bpozy);
  _leftEv(1,6) = ovsqrbar*(-bzopz + bpozz);

  _leftEv(2,1) = ov2a2*(-ascsnx - k2px - k2zx);
  _leftEv(2,2) = ov2a2*(-ascsny - k2py - k2zy);
  _leftEv(2,3) = ov2a2*(-ascsnz - k2pz - k2zz);
  _leftEv(2,4) = ov2sqra*(-afbppx - afbzzx);
  _leftEv(2,5) = ov2sqra*(-afbppy - afbzzy);
  _leftEv(2,6) = ov2sqra*(-afbppz - afbzzz);
  _leftEv(2,7) = ov2a2*alphas/rbar;

  _leftEv(3,0) = 1.;
  _leftEv(3,7) = -1./a2;

  _leftEv(4,4) = nx;
  _leftEv(4,5) = ny;
  _leftEv(4,6) = nz;

  _leftEv(5,1) = ov2a2*(ascsnx + k2px + k2zx);
  _leftEv(5,2) = ov2a2*(ascsny + k2py + k2zy);
  _leftEv(5,3) = ov2a2*(ascsnz + k2pz + k2zz);
  _leftEv(5,4) = ov2sqra*(-afbppx - afbzzx);
  _leftEv(5,5) = ov2sqra*(-afbppy - afbzzy);
  _leftEv(5,6) = ov2sqra*(-afbppz - afbzzz);
  _leftEv(5,7) = ov2a2*alphas/rbar;

  _leftEv(6,1) = -bzopx + bpozx;
  _leftEv(6,2) = -bzopy + bpozy;
  _leftEv(6,3) = -bzopz + bpozz;
  _leftEv(6,4) = ovsqrbar*(bzopx - bpozx);
  _leftEv(6,5) = ovsqrbar*(bzopy - bpozy);
  _leftEv(6,6) = ovsqrbar*(bzopz - bpozz);

  _leftEv(7,1) = ov2a2*(afcfnx - k1px - k1zx);
  _leftEv(7,2) = ov2a2*(afcfny - k1py - k1zy);
  _leftEv(7,3) = ov2a2*(afcfnz - k1pz - k1zz);
  _leftEv(7,4) = ov2sqra*(asbppx + asbzzx);
  _leftEv(7,5) = ov2sqra*(asbppy + asbzzy);
  _leftEv(7,6) = ov2sqra*(asbppz + asbzzz);
  _leftEv(7,7) = ov2a2*alphaf/rbar;

  //  CFLogDebugMax( "R*L = " << "\n" << RealMatrix(_rightEv*_leftEv) << "\n");
  // CFLogDebugMax( "L*R = " << "\n" << RealMatrix(_leftEv*_rightEv) << "\n");

  // compute the eigen values + and -
  for (CFuint iEq = 0; iEq < 8; ++iEq) {
    _eValuesP[iEq] = max((CFreal)0.,eValues[iEq]);
    _eValuesM[iEq] = min((CFreal)0.,eValues[iEq]);
  }

  // compute jacobian + and -
  jacobPlus = _rightEv*(_eValuesP*_leftEv);
  jacobMin  = _rightEv*(_eValuesM*_leftEv);
 }

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::computePhysicalData(const State& state,
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

  data[MHDTerm::VX] = u;
  data[MHDTerm::VY] = v;
  data[MHDTerm::VZ] = w;
  data[MHDTerm::BX] = Bx;
  data[MHDTerm::BY] = By;
  data[MHDTerm::BZ] = Bz;
  data[MHDTerm::RHO] = rho;
  data[MHDTerm::V] = sqrt(V2);
  data[MHDTerm::B] = sqrt(B2);
  data[MHDTerm::P] = state[7];
  data[MHDTerm::A] = sqrt(getModel()->getGamma()*data[MHDTerm::P]/rho);

  data[MHDTerm::GAMMA] = getModel()->getGamma();
  
  const RealVector& node = state.getCoordinates();
  data[MHDTerm::XP] =  node[XX];
  data[MHDTerm::YP] =  node[YY];
  data[MHDTerm::ZP] =  node[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::setDimensionalValuesPlusExtraValues(const State& state,
                                                              RealVector& result,
                                                              RealVector& extra)
{
  RealVector&  BDipole = getMagneticDipole(state.getCoordinates()[XX],
                                           state.getCoordinates()[YY],
                                           state.getCoordinates()[ZZ]);

  const CFreal gammaMinus1 = getModel()->getGamma() - 1.;

  const CFreal B1dotB0 = state[4]*BDipole[0] + state[5]*BDipole[1] + state[6]*BDipole[2];
  const CFreal sqB0 = BDipole[0]*BDipole[0] + BDipole[1]*BDipole[1] + BDipole[2]*BDipole[2];

  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal Bx1 = state[4];
  const CFreal By1 = state[5];
  const CFreal Bz1 = state[6];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal sqB1 = Bx1*Bx1 + By1*By1 + Bz1*Bz1;

  result[0] = state[0];
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
  result[4] = state[4];
  result[5] = state[5];
  result[6] = state[6];
  result[7] = state[7];

  extra.resize(9);

  if ((getModel()->getMX() != 0.0) || (getModel()->getMY() != 0.0) || (getModel()->getMZ() != 0.0)) {

     extra[0] = BDipole[0];
     extra[1] = BDipole[1];
     extra[2] = BDipole[2];
     extra[3] = state[4] + BDipole[0];
     extra[4] = state[5] + BDipole[1];
     extra[5] = state[6] + BDipole[2];
     extra[6] = sqrt(extra[3]*extra[3] + extra[4]*extra[4] + extra[5]*extra[5]);
     extra[7] = (state[7]/gammaMinus1) + 0.5*(rho*V2+sqB1) + B1dotB0 + 0.5*sqB0;

  }
  else {

      extra[0] = 0.0;
      extra[1] = 0.0;
      extra[2] = 0.0;
      extra[3] = state[4];
      extra[4] = state[5];
      extra[5] = state[6];
      extra[6] = sqrt(extra[3]*extra[3] + extra[4]*extra[4] + extra[5]*extra[5]);
      extra[7] = 0.0;

  }

  std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName = nsp + "_divBNodal";

  DataHandle<CFreal> divBNodal = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

  const CFuint stateID = state.getLocalID();

  extra[8] = divBNodal[stateID];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPrim::computeStateFromPhysicalData(const RealVector& data,
					State& state)
{
  const CFreal rho = data[MHDTerm::RHO];
  state[0] = rho;
  state[1] = data[MHDTerm::VX];
  state[2] = data[MHDTerm::VY];
  state[3] = data[MHDTerm::VZ];
  state[4] = data[MHDTerm::BX];
  state[5] = data[MHDTerm::BY];
  state[6] = data[MHDTerm::BZ];
  state[7] = data[MHDTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
