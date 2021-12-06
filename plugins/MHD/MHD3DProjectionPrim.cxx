#include "MHD/MHD.hh"

#include "MHD3DProjectionPrim.hh"
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

Environment::ObjectProvider<MHD3DProjectionPrim, ConvectiveVarSet, MHDModule, 1> 
mhd3DProjectionPrimProvider("MHD3DProjectionPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPrim::MHD3DProjectionPrim(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DProjectionVarSet(term),
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

MHD3DProjectionPrim::~MHD3DProjectionPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrim::setup()
{
  MHD3DProjectionVarSet::setup();

  _rm.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _rmInv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getNbEq());

}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD3DProjectionPrim::getExtraVarNames() const
{
  vector<std::string> names(9);
  names[0] = "BxPotential";
  names[1] = "ByPotential";
  names[2] = "BzPotential";
  names[3] = "BxTotal";
  names[4] = "ByTotal";
  names[5] = "BzTotal";
  names[6] = "BTotal";
  names[7] = "rhoETotal";
  //  names[8] = "divB";
  names[8] = "T";
  
  return names;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrim::computeJacobians()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrim::computeEigenValuesVectors(RealMatrix& rightEv,
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

  const CFreal rbar  = linearData[MHDProjectionTerm::RHO];
  const CFreal ubar  = linearData[MHDProjectionTerm::VX];
  const CFreal vbar  = linearData[MHDProjectionTerm::VY];
  const CFreal wbar  = linearData[MHDProjectionTerm::VZ];
  const CFreal Bxbar = linearData[MHDProjectionTerm::BX];
  const CFreal Bybar = linearData[MHDProjectionTerm::BY];
  const CFreal Bzbar = linearData[MHDProjectionTerm::BZ];
  const CFreal pbar  = linearData[MHDProjectionTerm::P];
  const CFreal phibar = linearData[MHDProjectionTerm::PHI];
  const CFreal Vn = ubar*nx  + vbar*ny + wbar*nz;
  // const CFreal Vp = ubar*px  + vbar*py + wbar*pz;
  // const CFreal Vz = ubar*zx  + vbar*zy + wbar*zz;

  const CFreal Bn = Bxbar*nx + Bybar*ny + Bzbar*nz;;
  const CFreal Bp = Bxbar*px + Bybar*py + Bzbar*pz;
  const CFreal Bz = Bxbar*zx + Bybar*zy + Bzbar*zz;

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
  // there was a bug in calculating the Alfven wave speed of sound that should be tested

  //const CFreal ca = Bn/sqrbar;
  const CFreal ca = std::abs(Bn)/sqrbar;

  CFreal a  = linearData[MHDProjectionTerm::A];
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

  if ((Bp*Bp + Bz*Bz) < MathTools::MathConsts::CFrealEps()) {
    bethap = 1.;
    bethaz = 0.;
  }
  else {
    bethap = Bp/sqrt(Bp*Bp + Bz*Bz);
    bethaz = Bz/sqrt(Bp*Bp + Bz*Bz);
  }

  const CFreal sgnbn = MathFunctions::sign(Bn);

  // sgnbn = 1. to make MHD work with B = 0

  const CFreal ovsq2 = 1./sqrt(2.);
  const CFreal ovsq2x = ovsq2*nx;
  const CFreal ovsq2y = ovsq2*ny;
  const CFreal ovsq2z = ovsq2*nz;

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
  // const CFreal afbz = alphaf*cf*bethaz*sgnbn;
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
  // const CFreal bzsqra = sqra*bethaz;
  // const CFreal ovbz = ov2sqra*bethaz;

  // matrix of right eigenvectors
  rightEv(0,0) = rbar*alphaf;
  rightEv(0,2) = rbar*alphas;
  rightEv(0,4) = 1.;
  rightEv(0,6) = rbar*alphas;
  rightEv(0,8) = rbar*alphaf;

  rightEv(1,0) = -afcfnx + k1px + k1zx;
  rightEv(1,1) = -bzopx + bethap*ovsq2*zx;
  rightEv(1,2) = -ascsnx - k2px - k2zx;
  rightEv(1,6) = ascsnx + k2px + k2zx;
  rightEv(1,7) = -bzopx + bethap*ovsq2*zx;
  rightEv(1,8) = afcfnx - k1px - k1zx;

  rightEv(2,0) = -afcfny + k1py + k1zy;
  rightEv(2,1) = -bzopy + bethap*ovsq2*zy;
  rightEv(2,2) = -ascsny - k2py - k2zy;
  rightEv(2,6) = ascsny + k2py + k2zy;
  rightEv(2,7) = -bzopy + bethap*ovsq2*zy;
  rightEv(2,8) = afcfny - k1py - k1zy;

  rightEv(3,0) = -afcfnz + k1pz + asbz*zz;
  rightEv(3,1) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  rightEv(3,2) = -alphas*cs*nz - k2pz - k2zz;
  rightEv(3,6) = alphas*cs*nz + k2pz + k2zz;
  rightEv(3,7) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  rightEv(3,8) = alphaf*cf*nz - k1pz - k1zz;

  rightEv(4,0) = asbppx*sqra + asbzzx*sqra;
  rightEv(4,1) = -sqrbar*bzopx + sqrbar*ovsq2*bethap*zx;
  rightEv(4,2) = -afbppx*sqra - afbzzx*sqra;
  rightEv(4,3) = ovsq2x;
  rightEv(4,5) = ovsq2x;
  rightEv(4,6) = -afbppx*sqra - afbzzx*sqra;
  rightEv(4,7) = sqrbar*bzopx - sqrbar*ovsq2*bethap*zx;
  rightEv(4,8) = asbppx*sqra + asbzzx*sqra;

  rightEv(5,0) = asbppy*sqra + asbzzy*sqra;
  rightEv(5,1) = -sqrbar*bzopy + sqrbar*ovsq2*bethap*zy;
  rightEv(5,2) = -afbppy*sqra - afbzzy*sqra;
  rightEv(5,3) = ovsq2y;
  rightEv(5,5) = ovsq2y;
  rightEv(5,6) = -afbppy*sqra - afbzzy*sqra;
  rightEv(5,7) = sqrbar*bzopy - sqrbar*ovsq2*bethap*zy;
  rightEv(5,8) = asbppy*sqra + asbzzy*sqra;

  rightEv(6,0) = asbppz*sqra + asbzzz*sqra;
  rightEv(6,1) = -sqrbar*bzopz + sqrbar*ovsq2*bethap*zz;
  rightEv(6,2) = -afbppz*sqra - afbzzz*sqra;
  rightEv(6,3) = ovsq2z;
  rightEv(6,5) = ovsq2z;
  rightEv(6,6) = -afbppz*sqra - afbzzz*sqra;
  rightEv(6,7) = sqrbar*bzopz - sqrbar*ovsq2*bethap*zz;
  rightEv(6,8) = asbppz*sqra + asbzzz*sqra;

  rightEv(7,0) = alphaf*gamma*pbar;
  rightEv(7,2) = alphas*gamma*pbar;
  rightEv(7,6) = alphas*gamma*pbar;
  rightEv(7,8) = alphaf*gamma*pbar;

  rightEv(8,3) = -refSpeed*ovsq2;
  rightEv(8,5) = refSpeed*ovsq2;

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

  leftEv(3,4) = ovsq2x;
  leftEv(3,5) = ovsq2y;
  leftEv(3,6) = ovsq2z;
  leftEv(3,8) = -1./refSpeed*ovsq2;

  leftEv(4,0) = 1.;
  leftEv(4,7) = -1./a2;

  leftEv(5,4) = ovsq2x;
  leftEv(5,5) = ovsq2y;
  leftEv(5,6) = ovsq2z;
  leftEv(5,8) = 1./refSpeed*ovsq2;

  leftEv(6,1) = ov2a2*(ascsnx + k2px + k2zx);
  leftEv(6,2) = ov2a2*(ascsny + k2py + k2zy);
  leftEv(6,3) = ov2a2*(ascsnz + k2pz + k2zz);
  leftEv(6,4) = ov2sqra*(-afbppx - afbzzx);
  leftEv(6,5) = ov2sqra*(-afbppy - afbzzy);
  leftEv(6,6) = ov2sqra*(-afbppz - afbzzz);
  leftEv(6,7) = ov2a2*alphas/rbar;

  leftEv(7,1) = -bzopx + bpozx;
  leftEv(7,2) = -bzopy + bpozy;
  leftEv(7,3) = -bzopz + bpozz;
  leftEv(7,4) = ovsqrbar*(bzopx - bpozx);
  leftEv(7,5) = ovsqrbar*(bzopy - bpozy);
  leftEv(7,6) = ovsqrbar*(bzopz - bpozz);

  leftEv(8,1) = ov2a2*(afcfnx - k1px - k1zx);
  leftEv(8,2) = ov2a2*(afcfny - k1py - k1zy);
  leftEv(8,3) = ov2a2*(afcfnz - k1pz - k1zz);
  leftEv(8,4) = ov2sqra*(asbppx + asbzzx);
  leftEv(8,5) = ov2sqra*(asbppy + asbzzy);
  leftEv(8,6) = ov2sqra*(asbppz + asbzzz);
  leftEv(8,7) = ov2a2*alphaf/rbar;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrim::setDimensionalValuesPlusExtraValues(const State& state,
                                                              RealVector& result,
                                                              RealVector& extra)
{
  const std::string potentialBType = getModel()->getPotentialBType();

  CFreal Bx0=0.0, By0=0.0, Bz0=0.0;
  RealVector stateCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
        stateCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
        BPFSSCartesian(PhysicalModelStack::getActive()->getDim());
  RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
        sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());

  if (potentialBType == "Dipole") {
     RealVector&  BDipole = getMagneticDipole(state.getCoordinates()[XX],
                                              state.getCoordinates()[YY],
                                              state.getCoordinates()[ZZ]);

     Bx0 = BDipole[0];
     By0 = BDipole[1];
     Bz0 = BDipole[2];
  }
  if (potentialBType == "PFSS") {
     stateCoordsCartesian[0] = state.getCoordinates()[XX];
     stateCoordsCartesian[1] = state.getCoordinates()[YY];
     stateCoordsCartesian[2] = state.getCoordinates()[ZZ];
     setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);
     computePFSSMagneticField(stateCoordsSpherical,BPFSSCartesian,sphCarTransMat);

     // Boltzmann constant
     const CFreal k = 1.3806503e-23;

     // magnetic permeability at vacuum
     const CFreal mu0 = 4.e-7*MathTools::MathConsts::CFrealPi();

     // mass of proton
     const CFreal mp = 1.67262158e-27;

     // mass of electron
     const CFreal me = 9.10938188e-31;

     const CFreal nRef = getModel()->getNRef();
     const CFreal TRef = getModel()->getTRef();

     const CFreal rhoRef = nRef*(mp+me);
     const CFreal vRef = sqrt(2.0*k*TRef/mp);

     const CFreal BRef = sqrt(mu0*rhoRef*vRef*vRef)*1.e4;

     Bx0 = BPFSSCartesian[0]/BRef;
     By0 = BPFSSCartesian[1]/BRef;
     Bz0 = BPFSSCartesian[2]/BRef;
  }

  const CFreal gammaMinus1 = getModel()->getGamma() - 1.;

  const CFreal B1dotB0 = state[4]*Bx0 + state[5]*By0 + state[6]*Bz0;
  const CFreal sqB0 = Bx0*Bx0 + By0*By0 + Bz0*Bz0;

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
  result[8] = state[8];

  extra.resize(9);

  extra[0] = Bx0;
  extra[1] = By0;
  extra[2] = Bz0;
  extra[3] = state[4] + extra[0];
  extra[4] = state[5] + extra[1];
  extra[5] = state[6] + extra[2];
  extra[6] = sqrt(extra[3]*extra[3] + extra[4]*extra[4] + extra[5]*extra[5]);
  extra[7] = (state[7]/gammaMinus1) + 0.5*(rho*V2+sqB1) + B1dotB0 + 0.5*sqB0;

  /*
  // AL: this code is buggy assuming nodal divB while data here can be cell-centered, causing indexing to fail!!!    
  std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName = nsp + "_divBNodal";
  DataHandle<CFreal> divBNodal = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);
  const CFuint stateID = state.getLocalID();
  extra[8] = divBNodal[stateID];
  */

  // calculate temperature here from perfect gas
  extra[8] = 0.0 ; // 
  
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrim::setRotationMatrices(const RealVector& normals)
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
  _rm(8,8) = 1.;


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
  _rmInv(8,8) = 1.;
}

//////////////////////////////////////////////////////////////////////////////

CFuint MHD3DProjectionPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrim::splitJacobian(RealMatrix& jacobPlus,
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

  const CFreal rbar  = linearData[MHDProjectionTerm::RHO];
  const CFreal ubar  = linearData[MHDProjectionTerm::VX];
  const CFreal vbar  = linearData[MHDProjectionTerm::VY];
  const CFreal wbar  = linearData[MHDProjectionTerm::VZ];
  const CFreal Bxbar = linearData[MHDProjectionTerm::BX];
  const CFreal Bybar = linearData[MHDProjectionTerm::BY];
  const CFreal Bzbar = linearData[MHDProjectionTerm::BZ];
  const CFreal pbar  = linearData[MHDProjectionTerm::P];
  const CFreal phibar = linearData[MHDProjectionTerm::PHI];
  const CFreal Vn = ubar*nx  + vbar*ny + wbar*nz;
  const CFreal Bn = Bxbar*nx + Bybar*ny + Bzbar*nz;;
  const CFreal Bp = Bxbar*px + Bybar*py + Bzbar*pz;
  const CFreal Bz = Bxbar*zx + Bybar*zy + Bzbar*zz;
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

  if ((Bp*Bp + Bz*Bz) < MathTools::MathConsts::CFrealEps()) {
    bethap = 1.;
    bethaz = 0.;
  }
  else {
    bethap = Bp/sqrt(Bp*Bp + Bz*Bz);
    bethaz = Bz/sqrt(Bp*Bp + Bz*Bz);
  }

  const CFreal sgnbn = MathFunctions::sign(Bn);

  // sgnbn = 1. to make MHD work with B = 0

  const CFreal ovsq2 = 1./sqrt(2.);
  const CFreal ovsq2x = ovsq2*nx;
  const CFreal ovsq2y = ovsq2*ny;
  const CFreal ovsq2z = ovsq2*nz;

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
  // const CFreal afbz = alphaf*cf*bethaz*sgnbn;
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
  // const CFreal bzsqra = sqra*bethaz;
  // const CFreal ovbz = ov2sqra*bethaz;

  // matrix of right eigenvectors
  _rightEv(0,0) = rbar*alphaf;
  _rightEv(0,2) = rbar*alphas;
  _rightEv(0,4) = 1.;
  _rightEv(0,6) = rbar*alphas;
  _rightEv(0,8) = rbar*alphaf;

  _rightEv(1,0) = -afcfnx + k1px + k1zx;
  _rightEv(1,1) = -bzopx + bethap*ovsq2*zx;
  _rightEv(1,2) = -ascsnx - k2px - k2zx;
  _rightEv(1,6) = ascsnx + k2px + k2zx;
  _rightEv(1,7) = -bzopx + bethap*ovsq2*zx;
  _rightEv(1,8) = afcfnx - k1px - k1zx;

  _rightEv(2,0) = -afcfny + k1py + k1zy;
  _rightEv(2,1) = -bzopy + bethap*ovsq2*zy;
  _rightEv(2,2) = -ascsny - k2py - k2zy;
  _rightEv(2,6) = ascsny + k2py + k2zy;
  _rightEv(2,7) = -bzopy + bethap*ovsq2*zy;
  _rightEv(2,8) = afcfny - k1py - k1zy;

  _rightEv(3,0) = -afcfnz + k1pz + asbz*zz;
  _rightEv(3,1) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  _rightEv(3,2) = -alphas*cs*nz - k2pz - k2zz;
  _rightEv(3,6) = alphas*cs*nz + k2pz + k2zz;
  _rightEv(3,7) = -bethaz*ovsq2*pz + bethap*ovsq2*zz;
  _rightEv(3,8) = alphaf*cf*nz - k1pz - k1zz;

  _rightEv(4,0) = asbppx*sqra + asbzzx*sqra;
  _rightEv(4,1) = -sqrbar*bzopx + sqrbar*ovsq2*bethap*zx;
  _rightEv(4,2) = -afbppx*sqra - afbzzx*sqra;
  _rightEv(4,3) = ovsq2x;
  _rightEv(4,5) = ovsq2x;
  _rightEv(4,6) = -afbppx*sqra - afbzzx*sqra;
  _rightEv(4,7) = sqrbar*bzopx - sqrbar*ovsq2*bethap*zx;
  _rightEv(4,8) = asbppx*sqra + asbzzx*sqra;

  _rightEv(5,0) = asbppy*sqra + asbzzy*sqra;
  _rightEv(5,1) = -sqrbar*bzopy + sqrbar*ovsq2*bethap*zy;
  _rightEv(5,2) = -afbppy*sqra - afbzzy*sqra;
  _rightEv(5,3) = ovsq2y;
  _rightEv(5,5) = ovsq2y;
  _rightEv(5,6) = -afbppy*sqra - afbzzy*sqra;
  _rightEv(5,7) = sqrbar*bzopy - sqrbar*ovsq2*bethap*zy;
  _rightEv(5,8) = asbppy*sqra + asbzzy*sqra;

  _rightEv(6,0) = asbppz*sqra + asbzzz*sqra;
  _rightEv(6,1) = -sqrbar*bzopz + sqrbar*ovsq2*bethap*zz;
  _rightEv(6,2) = -afbppz*sqra - afbzzz*sqra;
  _rightEv(6,3) = ovsq2z;
  _rightEv(6,5) = ovsq2z;
  _rightEv(6,6) = -afbppz*sqra - afbzzz*sqra;
  _rightEv(6,7) = sqrbar*bzopz - sqrbar*ovsq2*bethap*zz;
  _rightEv(6,8) = asbppz*sqra + asbzzz*sqra;

  _rightEv(7,0) = alphaf*gamma*pbar;
  _rightEv(7,2) = alphas*gamma*pbar;
  _rightEv(7,6) = alphas*gamma*pbar;
  _rightEv(7,8) = alphaf*gamma*pbar;

  _rightEv(8,3) = -refSpeed*ovsq2;
  _rightEv(8,5) = refSpeed*ovsq2;

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

  _leftEv(3,4) = ovsq2x;
  _leftEv(3,5) = ovsq2y;
  _leftEv(3,6) = ovsq2z;
  _leftEv(3,8) = -1./refSpeed*ovsq2;

  _leftEv(4,0) = 1.;
  _leftEv(4,7) = -1./a2;

  _leftEv(5,4) = ovsq2x;
  _leftEv(5,5) = ovsq2y;
  _leftEv(5,6) = ovsq2z;
  _leftEv(5,8) = 1./refSpeed*ovsq2;

  _leftEv(6,1) = ov2a2*(ascsnx + k2px + k2zx);
  _leftEv(6,2) = ov2a2*(ascsny + k2py + k2zy);
  _leftEv(6,3) = ov2a2*(ascsnz + k2pz + k2zz);
  _leftEv(6,4) = ov2sqra*(-afbppx - afbzzx);
  _leftEv(6,5) = ov2sqra*(-afbppy - afbzzy);
  _leftEv(6,6) = ov2sqra*(-afbppz - afbzzz);
  _leftEv(6,7) = ov2a2*alphas/rbar;

  _leftEv(7,1) = -bzopx + bpozx;
  _leftEv(7,2) = -bzopy + bpozy;
  _leftEv(7,3) = -bzopz + bpozz;
  _leftEv(7,4) = ovsqrbar*(bzopx - bpozx);
  _leftEv(7,5) = ovsqrbar*(bzopy - bpozy);
  _leftEv(7,6) = ovsqrbar*(bzopz - bpozz);

  _leftEv(8,1) = ov2a2*(afcfnx - k1px - k1zx);
  _leftEv(8,2) = ov2a2*(afcfny - k1py - k1zy);
  _leftEv(8,3) = ov2a2*(afcfnz - k1pz - k1zz);
  _leftEv(8,4) = ov2sqra*(asbppx + asbzzx);
  _leftEv(8,5) = ov2sqra*(asbppy + asbzzy);
  _leftEv(8,6) = ov2sqra*(asbppz + asbzzz);
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

void MHD3DProjectionPrim::computePhysicalData(const State& state,
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
  data[MHDTerm::ZP] =  node[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPrim::computeStateFromPhysicalData(const RealVector& data,
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

