#include "MHD/MHD.hh"
#include "MHD3DCons.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodData.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/NotImplementedException.hh"
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

Environment::ObjectProvider<MHD3DCons, ConvectiveVarSet, MHDModule, 1> mhd3DConsProvider("MHD3DCons");

//////////////////////////////////////////////////////////////////////////////

MHD3DCons::MHD3DCons(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DVarSet(term),
  _rm(),
  _rmInv(),
  _linearState(),
  _linState(),
  _stateFlux(),
  _dudw(),
  _rightMat(),
  _leftStateVector(),
  _rightStateVector(),
  _leftStateVectorRot(),
  _rightStateVectorRot(),
  _dw()
{
  vector<std::string> names(8);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoW";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "rhoE";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DCons::~MHD3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::setup()
{
  MHD3DVarSet::setup();

  _rm.resize(PhysicalModelStack::getActive()->getNbEq(),
      PhysicalModelStack::getActive()->getNbEq());
  _rmInv.resize(PhysicalModelStack::getActive()->getNbEq(),
   PhysicalModelStack::getActive()->getNbEq());
  _linearState.resize(PhysicalModelStack::getActive()->getNbEq());
  _linState.resize(PhysicalModelStack::getActive()->getNbEq());
  _stateFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _dudw.resize(PhysicalModelStack::getActive()->getNbEq(),
  PhysicalModelStack::getActive()->getNbEq());
  _rightMat.resize(PhysicalModelStack::getActive()->getNbEq(),
      PhysicalModelStack::getActive()->getNbEq());

  _leftStateVector.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _rightStateVector.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _leftStateVectorRot.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _rightStateVectorRot.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
  _dw.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD3DCons::getExtraVarNames() const
{
  vector<std::string> names(10);
  names[0] = "BxDipole";
  names[1] = "ByDipole";
  names[2] = "BzDipole";
  names[3] = "BxTotal";
  names[4] = "ByTotal";
  names[5] = "BzTotal";
  names[6] = "BTotal";
  names[7] = "rhoETotal";
  names[8] = "p";
  names[9] = "divB";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"MHD3DCons::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::computeEigenValuesVectors(RealMatrix& rightEv,
              RealMatrix& leftEv,
              RealVector& eValues,
              const RealVector& normal)
{
  // components of the rotated linearized state vector
  setRotationMatrices(normal, _rm, _rmInv);

  const RealVector& linearData = getModel()->getPhysicalData();

  //set the linear states in cons variables
  const CFreal rho = linearData[MHDTerm::RHO];
  _linearState[0] = rho;
  _linearState[1] = rho*linearData[MHDTerm::VX];
  _linearState[2] = rho*linearData[MHDTerm::VY];
  _linearState[3] = rho*linearData[MHDTerm::VZ];
  _linearState[4] = linearData[MHDTerm::BX];
  _linearState[5] = linearData[MHDTerm::BY];
  _linearState[6] = linearData[MHDTerm::BZ];
  const CFreal magB = linearData[MHDTerm::B];
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.);
  const CFreal rhoH = gammaDivGammaMinus1*linearData[MHDTerm::P]
    + 0.5*rho*linearData[MHDTerm::V]*linearData[MHDTerm::V] + 0.5*magB*magB;
  _linearState[7] = rhoH - linearData[MHDTerm::P];

  _linState = _rm*_linearState;

  const CFreal gammaMinus1 = getModel()->getGamma() - 1.0;
  const CFreal rhoAv  = _linState[0];
  const CFreal vnAv  = _linState[1]/rhoAv;
  const CFreal vtAv  = _linState[2]/rhoAv;
  const CFreal vzAv  = _linState[3]/rhoAv;
  const CFreal BnAv = _linState[4];
  const CFreal BtAv = _linState[5];
  const CFreal BzAv = _linState[6];
  const CFreal v2Av = vnAv*vnAv + vtAv*vtAv + vzAv*vzAv;
  const CFreal B2Av = BnAv*BnAv + BtAv*BtAv + BzAv*BzAv;
  const CFreal pAv  = gammaMinus1*(_linState[7] - 0.5*rhoAv*v2Av
       - 0.5*B2Av);

  // speeds of sound of the system

  const CFreal caAv = std::abs(BnAv)/sqrt(rhoAv);

  const CFreal astar2 = (gamma*pAv + B2Av)/rhoAv;
  CFreal cfAv2 = 0.5*(astar2 + sqrt(astar2*astar2
              - 4.0*gamma*pAv*BnAv*BnAv/(rhoAv*rhoAv)));
  CFreal csAv2 = 0.5*(astar2 - sqrt(astar2*astar2
                - 4.0*gamma*pAv*BnAv*BnAv/(rhoAv*rhoAv)));

  //if (std::abs(cfAv2) < 1.0e-16) cfAv2 = 0.;

  // absolute value is necessary to avoid from negative csAv2

  const CFreal cfAv = sqrt(cfAv2);
  const CFreal csAv = sqrt(std::abs(csAv2));

  const CFreal aAv  = sqrt(gamma*pAv/rhoAv);
  CFreal aAv2 = aAv*aAv;

  // eigenvalues of the system

  eValues[0] = vnAv;
  eValues[1] = vnAv;
  eValues[2] = vnAv + caAv;
  eValues[3] = vnAv - caAv;
  eValues[4] = vnAv + cfAv;
  eValues[5] = vnAv - cfAv;
  eValues[6] = vnAv + csAv;
  eValues[7] = vnAv - csAv;

  // Scaling factors for the normalized eigensystem and the
  // singularities in eigenvectors

  // fix against machine precision
 //  if (csAv2 > aAv2) {
//     aAv2 = csAv2;
//   }
//   if (cfAv2 < aAv2) {
//     aAv2 = cfAv2;
//   }

//   CFreal alphaf2 = (aAv2 - csAv2)/(cfAv2 - csAv2);
//   CFreal alphas2 = (cfAv2 - aAv2)/(cfAv2 - csAv2);
//   if (std::abs(alphaf2) < MathTools::MathConsts::CFrealEps()) alphaf2 = 0.;
//   if (std::abs(alphas2) < MathTools::MathConsts::CFrealEps()) alphas2 = 0.;


//   const CFreal alphaf = sqrt(alphaf2);
//   const CFreal alphas = sqrt(alphas2);

  CFreal alphaf = 0.0;
  CFreal alphas = 0.0;

  if (std::abs(cfAv2-csAv2) < MathTools::MathConsts::CFrealEps())
    if (sqrt(BtAv*BtAv + BzAv*BzAv) == 0.0) {

      alphaf = MathTools::MathFunctions::heavyside(aAv-caAv);
      alphas = MathTools::MathFunctions::heavyside(caAv-aAv);

    }
    else {

      const CFreal phi = atan((caAv-aAv)/sqrt(BtAv*BtAv + BzAv*BzAv));
      alphaf = sin(phi/2.0);
      alphas = cos(phi/2.0);

    }
  else {

    alphaf = sqrt(std::abs((aAv2 - csAv2)/(cfAv2 - csAv2)));
    alphas = sqrt(std::abs((cfAv2 - aAv2)/(cfAv2 - csAv2)));

  }

  CFreal bethap = 0.0;
  CFreal bethaz = 0.0;

  const CFreal sgnbn = MathTools::MathFunctions::signum(BnAv);
  const CFreal ovsq2 = 1.0/sqrt(2.0);

  if (sqrt(BtAv*BtAv + BzAv*BzAv) < MathTools::MathConsts::CFrealEps()) {

    bethap = ovsq2;
    bethaz = ovsq2;

  }
  else {

    bethap = BtAv/sqrt(BtAv*BtAv + BzAv*BzAv);
    bethaz = BzAv/sqrt(BtAv*BtAv + BzAv*BzAv);

  }

  // right eigenvector matrix in primitive variables

  _rightMat(0,0) = 1.0;
  _rightMat(1,0) = 0.0;
  _rightMat(2,0) = 0.0;
  _rightMat(3,0) = 0.0;
  _rightMat(4,0) = 0.0;
  _rightMat(5,0) = 0.0;
  _rightMat(6,0) = 0.0;
  _rightMat(7,0) = 0.0;

  _rightMat(0,1) = 0.0;
  _rightMat(1,1) = 0.0;
  _rightMat(2,1) = 0.0;
  _rightMat(3,1) = 0.0;
  _rightMat(4,1) = 1.0;
  _rightMat(5,1) = 0.0;
  _rightMat(6,1) = 0.0;
  _rightMat(7,1) = 0.0;

  _rightMat(0,2) = 0.0;
  _rightMat(1,2) = 0.0;
  _rightMat(2,2) = -bethaz*ovsq2;
  _rightMat(3,2) = bethap*ovsq2;
  _rightMat(4,2) = 0.0;
  _rightMat(5,2) = sqrt(rhoAv)*bethaz*ovsq2;
  _rightMat(6,2) = -sqrt(rhoAv)*bethap*ovsq2;
  _rightMat(7,2) = 0.0;

  _rightMat(0,3) = 0.0;
  _rightMat(1,3) = 0.0;
  _rightMat(2,3) = -bethaz*ovsq2;
  _rightMat(3,3) = bethap*ovsq2;
  _rightMat(4,3) = 0.0;
  _rightMat(5,3) = -sqrt(rhoAv)*bethaz*ovsq2;
  _rightMat(6,3) = sqrt(rhoAv)*bethap*ovsq2;
  _rightMat(7,3) = 0.0;

  _rightMat(0,4) = rhoAv*alphaf;
  _rightMat(1,4) = alphaf*cfAv;
  _rightMat(2,4) = -alphas*csAv*bethap*sgnbn;
  _rightMat(3,4) = -alphas*csAv*bethaz*sgnbn;
  _rightMat(4,4) = 0.0;
  _rightMat(5,4) = alphas*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,4) = alphas*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,4) = alphaf*gamma*pAv;

  _rightMat(0,5) = rhoAv*alphaf;
  _rightMat(1,5) = -alphaf*cfAv;
  _rightMat(2,5) = alphas*csAv*bethap*sgnbn;
  _rightMat(3,5) = alphas*csAv*bethaz*sgnbn;
  _rightMat(4,5) = 0.0;
  _rightMat(5,5) = alphas*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,5) = alphas*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,5) = alphaf*gamma*pAv;

  _rightMat(0,6) = rhoAv*alphas;
  _rightMat(1,6) = alphas*csAv;
  _rightMat(2,6) = alphaf*cfAv*bethap*sgnbn;
  _rightMat(3,6) = alphaf*cfAv*bethaz*sgnbn;
  _rightMat(4,6) = 0.0;
  _rightMat(5,6) = -alphaf*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,6) = -alphaf*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,6) = alphas*gamma*pAv;

  _rightMat(0,7) = rhoAv*alphas;
  _rightMat(1,7) = -alphas*csAv;
  _rightMat(2,7) = -alphaf*cfAv*bethap*sgnbn;
  _rightMat(3,7) = -alphaf*cfAv*bethaz*sgnbn;
  _rightMat(4,7) = 0.0;
  _rightMat(5,7) = -alphaf*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,7) = -alphaf*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,7) = alphas*gamma*pAv;

  // the matrix of dU/dW where U is the conservative variables and
  // W is the primitive variables

  _dudw(0,0) = 1.0;
  _dudw(1,0) = vnAv;
  _dudw(2,0) = vtAv;
  _dudw(3,0) = vzAv;
  _dudw(4,0) = 0.0;
  _dudw(5,0) = 0.0;
  _dudw(6,0) = 0.0;
  _dudw(7,0) = 0.5*v2Av;

  _dudw(0,1) = 0.0;
  _dudw(1,1) = rhoAv;
  _dudw(2,1) = 0.0;
  _dudw(3,1) = 0.0;
  _dudw(4,1) = 0.0;
  _dudw(5,1) = 0.0;
  _dudw(6,1) = 0.0;
  _dudw(7,1) = rhoAv*vnAv;

  _dudw(0,2) = 0.0;
  _dudw(1,2) = 0.0;
  _dudw(2,2) = rhoAv;
  _dudw(3,2) = 0.0;
  _dudw(4,2) = 0.0;
  _dudw(5,2) = 0.0;
  _dudw(6,2) = 0.0;
  _dudw(7,2) = rhoAv*vtAv;

  _dudw(0,3) = 0.0;
  _dudw(1,3) = 0.0;
  _dudw(2,3) = 0.0;
  _dudw(3,3) = rhoAv;
  _dudw(4,3) = 0.0;
  _dudw(5,3) = 0.0;
  _dudw(6,3) = 0.0;
  _dudw(7,3) = rhoAv*vzAv;

  _dudw(0,4) = 0.0;
  _dudw(1,4) = 0.0;
  _dudw(2,4) = 0.0;
  _dudw(3,4) = 0.0;
  _dudw(4,4) = 1.0;
  _dudw(5,4) = 0.0;
  _dudw(6,4) = 0.0;
  _dudw(7,4) = BnAv;

  _dudw(0,5) = 0.0;
  _dudw(1,5) = 0.0;
  _dudw(2,5) = 0.0;
  _dudw(3,5) = 0.0;
  _dudw(4,5) = 0.0;
  _dudw(5,5) = 1.0;
  _dudw(6,5) = 0.0;
  _dudw(7,5) = BtAv;

  _dudw(0,6) = 0.0;
  _dudw(1,6) = 0.0;
  _dudw(2,6) = 0.0;
  _dudw(3,6) = 0.0;
  _dudw(4,6) = 0.0;
  _dudw(5,6) = 0.0;
  _dudw(6,6) = 1.0;
  _dudw(7,6) = BzAv;

  _dudw(0,7) = 0.0;
  _dudw(1,7) = 0.0;
  _dudw(2,7) = 0.0;
  _dudw(3,7) = 0.0;
  _dudw(4,7) = 0.0;
  _dudw(5,7) = 0.0;
  _dudw(6,7) = 0.0;
  _dudw(7,7) = 1.0/(gamma-1.0);

  // right eigenvector matrix in conservative variables

  rightEv = _dudw*_rightMat;

  // left eigenvector matrix in primitive variables

  leftEv(0,0) = 1.0;
  leftEv(0,1) = 0.0;
  leftEv(0,2) = 0.0;
  leftEv(0,3) = 0.0;
  leftEv(0,4) = 0.0;
  leftEv(0,5) = 0.0;
  leftEv(0,6) = 0.0;
  leftEv(0,7) = -1.0/aAv2;

  leftEv(1,0) = 0.0;
  leftEv(1,1) = 0.0;
  leftEv(1,2) = 0.0;
  leftEv(1,3) = 0.0;
  leftEv(1,4) = 1.0;
  leftEv(1,5) = 0.0;
  leftEv(1,6) = 0.0;
  leftEv(1,7) = 0.0;

  leftEv(2,0) = 0.0;
  leftEv(2,1) = 0.0;
  leftEv(2,2) = -bethaz*ovsq2;
  leftEv(2,3) = bethap*ovsq2;
  leftEv(2,4) = 0.0;
  leftEv(2,5) = bethaz*ovsq2/sqrt(rhoAv);
  leftEv(2,6) = -bethap*ovsq2/sqrt(rhoAv);
  leftEv(2,7) = 0.0;

  leftEv(3,0) = 0.0;
  leftEv(3,1) = 0.0;
  leftEv(3,2) = -bethaz*ovsq2;
  leftEv(3,3) = bethap*ovsq2;
  leftEv(3,4) = 0.0;
  leftEv(3,5) = -bethaz*ovsq2/sqrt(rhoAv);
  leftEv(3,6) = bethap*ovsq2/sqrt(rhoAv);
  leftEv(3,7) = 0.0;

  leftEv(4,0) = 0.0;
  leftEv(4,1) = alphaf*cfAv/(2.0*aAv2);
  leftEv(4,2) = -alphas*csAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(4,3) = -alphas*csAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(4,4) = 0.0;
  leftEv(4,5) = alphas*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(4,6) = alphas*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(4,7) = alphaf/(2.0*rhoAv*aAv2);

  leftEv(5,0) = 0.0;
  leftEv(5,1) = -alphaf*cfAv/(2.0*aAv2);
  leftEv(5,2) = alphas*csAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(5,3) = alphas*csAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(5,4) = 0.0;
  leftEv(5,5) = alphas*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(5,6) = alphas*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(5,7) = alphaf/(2.0*rhoAv*aAv2);

  leftEv(6,0) = 0.0;
  leftEv(6,1) = alphas*csAv/(2.0*aAv2);
  leftEv(6,2) = alphaf*cfAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(6,3) = alphaf*cfAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(6,4) = 0.0;
  leftEv(6,5) = -alphaf*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(6,6) = -alphaf*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(6,7) = alphas/(2.0*rhoAv*aAv2);

  leftEv(7,0) = 0.0;
  leftEv(7,1) = -alphas*csAv/(2.0*aAv2);
  leftEv(7,2) = -alphaf*cfAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(7,3) = -alphaf*cfAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(7,4) = 0.0;
  leftEv(7,5) = -alphaf*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(7,6) = -alphaf*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(7,7) = alphas/(2.0*rhoAv*aAv2);
}

//////////////////////////////////////////////////////////////////////////////

CFuint MHD3DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::computePhysicalData(const State& state,
				    RealVector& data)
{
  const CFreal rho = state[0];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal w = state[3]/rho;
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;
  const CFreal gamma = getModel()->getGamma();

  data[MHDTerm::VX] = u;
  data[MHDTerm::VY] = v;
  data[MHDTerm::VZ] = w;
  data[MHDTerm::BX] = Bx;
  data[MHDTerm::BY] = By;
  data[MHDTerm::BZ] = Bz;
  data[MHDTerm::B] = sqrt(B2);
  data[MHDTerm::RHO] = rho;
  data[MHDTerm::V] = sqrt(V2);
  data[MHDTerm::P] = (gamma - 1.)*(state[7] - 0.5*(rho*V2 + B2));

  // this cures problems due to occasionally negative pressures
  if (data[MHDTerm::P] < 0.) {
    CFLog(VERBOSE, "ERROR: p<0 (= " << data[MHDTerm::P] << ")\n");
    data[MHDTerm::P] = 0.;
  }

  const CFreal a2 = gamma*data[MHDTerm::P]/rho;
  data[MHDTerm::A] = sqrt(a2);

  data[MHDTerm::GAMMA] = gamma;
 
  const RealVector& node = state.getCoordinates();
  data[MHDTerm::XP] =  node[XX];
  data[MHDTerm::YP] =  node[YY];
  data[MHDTerm::ZP] =  node[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::setRotationMatrices(const RealVector& normals,
            RealMatrix& rm,
            RealMatrix& rmInv)
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

  rm(0,0) = 1.;
  rm(1,1) = nx;
  rm(1,2) = ny;
  rm(1,3) = nz;

  rm(2,1) = px;
  rm(2,2) = py;
  rm(2,3) = pz;

  rm(3,1) = zx;
  rm(3,2) = zy;
  rm(3,3) = zz;

  rm(4,4) = nx;
  rm(4,5) = ny;
  rm(4,6) = nz;

  rm(5,4) = px;
  rm(5,5) = py;
  rm(5,6) = pz;

  rm(6,4) = zx;
  rm(6,5) = zy;
  rm(6,6) = zz;

  rm(7,7) = 1.;


  rmInv(0,0) = 1.;
  rmInv(1,1) = nx;
  rmInv(1,2) = px;
  rmInv(1,3) = zx;

  rmInv(2,1) = ny;
  rmInv(2,2) = py;
  rmInv(2,3) = zy;

  rmInv(3,1) = nz;
  rmInv(3,2) = pz;
  rmInv(3,3) = zz;

  rmInv(4,4) = nx;
  rmInv(4,5) = px;
  rmInv(4,6) = zx;

  rmInv(5,4) = ny;
  rmInv(5,5) = py;
  rmInv(5,6) = zy;

  rmInv(6,4) = nz;
  rmInv(6,5) = pz;
  rmInv(6,6) = zz;

  rmInv(7,7) = 1.;

  // // check for rm * rmInv = I

//   RealMatrix matrixMultip(PhysicalModelStack::getActive()->getNbEq(),
//        PhysicalModelStack::getActive()->getNbEq());

//   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

//   matrixMultip = rm*rmInv;

//   for(int i = 0; i < nbEqs; i++)
//     for(int j = 0; j < nbEqs; j++)
//       if (matrixMultip(i,j) != 1.0)
//  cout << "here " << matrixMultip(i,j) << endl;

  CFLogDebugMax( "RM*RMI = " << "\n" << RealMatrix(rm*rmInv) << "\n");
  CFLogDebugMax( "RMI*RM = " << "\n" << RealMatrix(rmInv*rm) << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::setWaveStrengths(RealVector& waveStrengths,
         const RealMatrix& leftEv,
         State& leftState,
         State& rightState,
         const RealVector& normal)
{
  // components of the rotated linearized state vector

  setRotationMatrices(normal, _rm, _rmInv);

  // leftState and rightState should be in rotated form
  // in order to be used for the variable _sumFlux
  // leftState and rightState are modified accordingly
  RealVector& leftStateVector = _leftStateVector;
  RealVector& rightStateVector = _rightStateVector;
  RealVector& leftStateVectorRot = _leftStateVectorRot;
  RealVector& rightStateVectorRot = _rightStateVectorRot;

  leftStateVector[0] = leftState[0];
  leftStateVector[1] = leftState[1];
  leftStateVector[2] = leftState[2];
  leftStateVector[3] = leftState[3];
  leftStateVector[4] = leftState[4];
  leftStateVector[5] = leftState[5];
  leftStateVector[6] = leftState[6];
  leftStateVector[7] = leftState[7];

  rightStateVector[0] = rightState[0];
  rightStateVector[1] = rightState[1];
  rightStateVector[2] = rightState[2];
  rightStateVector[3] = rightState[3];
  rightStateVector[4] = rightState[4];
  rightStateVector[5] = rightState[5];
  rightStateVector[6] = rightState[6];
  rightStateVector[7] = rightState[7];

  leftStateVectorRot = _rm*leftStateVector;
  rightStateVectorRot = _rm*rightStateVector;

  leftState[0] = leftStateVectorRot[0];
  leftState[1] = leftStateVectorRot[1];
  leftState[2] = leftStateVectorRot[2];
  leftState[3] = leftStateVectorRot[3];
  leftState[4] = leftStateVectorRot[4];
  leftState[5] = leftStateVectorRot[5];
  leftState[6] = leftStateVectorRot[6];
  leftState[7] = leftStateVectorRot[7];

  rightState[0] = rightStateVectorRot[0];
  rightState[1] = rightStateVectorRot[1];
  rightState[2] = rightStateVectorRot[2];
  rightState[3] = rightStateVectorRot[3];
  rightState[4] = rightStateVectorRot[4];
  rightState[5] = rightStateVectorRot[5];
  rightState[6] = rightStateVectorRot[6];
  rightState[7] = rightStateVectorRot[7];

  const CFreal gammaMinus1 = getModel()->getGamma() - 1.0;

  const CFreal rhoL = leftState[0];
  const CFreal rhoR = rightState[0];
  const CFreal vnL = leftState[1]/rhoL;
  const CFreal vnR = rightState[1]/rhoR;
  const CFreal vtL = leftState[2]/rhoL;
  const CFreal vtR = rightState[2]/rhoR;
  const CFreal vzL = leftState[3]/rhoL;
  const CFreal vzR = rightState[3]/rhoR;
  const CFreal BnL = leftState[4];
  const CFreal BnR = rightState[4];
  const CFreal BtL = leftState[5];
  const CFreal BtR = rightState[5];
  const CFreal BzL = leftState[6];
  const CFreal BzR = rightState[6];
  const CFreal rhoEL = leftState[7];
  const CFreal rhoER = rightState[7];
  const CFreal v2L = vnL*vnL + vtL*vtL + vzL*vzL;
  const CFreal v2R = vnR*vnR + vtR*vtR + vzR*vzR;
  const CFreal B2L = BnL*BnL + BtL*BtL + BzL*BzL;
  const CFreal B2R = BnR*BnR + BtR*BtR + BzR*BzR;

  const CFreal pL = gammaMinus1*(rhoEL - 0.5*rhoL*v2L - 0.5*B2L);
  const CFreal pR = gammaMinus1*(rhoER - 0.5*rhoR*v2R - 0.5*B2R);

  RealVector& dw = _dw;

  dw[0] = rhoR - rhoL;
  dw[1] = vnR - vnL;
  dw[2] = vtR - vtL;
  dw[3] = vzR - vzL;
  dw[4] = BnR - BnL;
  dw[5] = BtR - BtL;
  dw[6] = BzR - BzL;
  dw[7] = pR - pL;

  waveStrengths = leftEv*dw;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::setDimensionalValuesPlusExtraValues(const State& state,
                RealVector& result,
                RealVector& extra)
{
  RealVector&  BDipole = getMagneticDipole(state.getCoordinates()[XX],
					   state.getCoordinates()[YY],
					   state.getCoordinates()[ZZ]);
  const CFreal B1dotB0 = state[4]*BDipole[0] + state[5]*BDipole[1] + state[6]*BDipole[2];
  const CFreal sqB0 = BDipole[0]*BDipole[0] + BDipole[1]*BDipole[1] + BDipole[2]*BDipole[2];
  const CFreal rho = state[0];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal w = state[3]/rho;
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal sqB1 = state[4]*state[4] + state[5]*state[5] + state[6]*state[6];

  result[0] = state[0];
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
  result[4] = state[4];
  result[5] = state[5];
  result[6] = state[6];
  result[7] = state[7];

  extra.resize(10);

  if ((getModel()->getMX() != 0.0) || (getModel()->getMY() != 0.0) || (getModel()->getMZ() != 0.0)) {
     extra[0] = BDipole[0];
     extra[1] = BDipole[1];
     extra[2] = BDipole[2];
     extra[3] = state[4] + BDipole[0];
     extra[4] = state[5] + BDipole[1];
     extra[5] = state[6] + BDipole[2];
     extra[6] = sqrt(extra[3]*extra[3] + extra[4]*extra[4] + extra[5]*extra[5]);
     extra[7] = state[7] + B1dotB0 + 0.5*sqB0;

  }

  else {

     extra[0] = 0.0;
     extra[1] = 0.0;
     extra[2] = 0.0;
     extra[3] = 0.0;
     extra[4] = 0.0;
     extra[5] = 0.0;
     extra[6] = 0.0;
     extra[7] = 0.0;

  }

  extra[8] = (getModel()->getGamma() - 1.)*
    (state[7] - 0.5*(rho*V2 + sqB1));

  std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName = nsp + "_divBNodal";

  DataHandle<CFreal> divBNodal = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

  const CFuint stateID = state.getLocalID();

  extra[9] = divBNodal[stateID];

}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHD3DCons::getRotatedFlux(const Framework::State& state,
              const RealVector& normal)
{
  const CFreal gammaMinus1 = getModel()->getGamma() - 1.0;

  const CFreal rhoL = state[0];
  const CFreal vnL = state[1]/rhoL;
  const CFreal vtL = state[2]/rhoL;
  const CFreal vzL = state[3]/rhoL;
  const CFreal BnL = state[4];
  const CFreal BtL = state[5];
  const CFreal BzL = state[6];
  const CFreal rhoEL = state[7];
  const CFreal v2L = vnL*vnL + vtL*vtL + vzL*vzL;
  const CFreal B2L = BnL*BnL + BtL*BtL + BzL*BzL;
  const CFreal pL = gammaMinus1*(rhoEL - 0.5*rhoL*v2L - 0.5*B2L);

  _stateFlux[0] = rhoL*vnL;
  _stateFlux[1] = rhoL*vnL*vnL + (pL+0.5*B2L) - BnL*BnL;
  _stateFlux[2] = rhoL*vnL*vtL - BnL*BtL;
  _stateFlux[3] = rhoL*vnL*vzL - BnL*BzL;
  _stateFlux[4] = 0.0;
  _stateFlux[5] = vnL*BtL - vtL*BnL;
  _stateFlux[6] = vnL*BzL - vzL*BnL;
  _stateFlux[7] = (rhoEL+pL+0.5*B2L)*vnL-BnL*(vnL*BnL+vtL*BtL+vzL*BzL);

  return _stateFlux;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCons::computeStateFromPhysicalData(const RealVector& data,
            State& state)
{
  const CFreal rho = data[MHDTerm::RHO];
  state[0] = rho;
  state[1] = rho*data[MHDTerm::VX];
  state[2] = rho*data[MHDTerm::VY];
  state[3] = rho*data[MHDTerm::VZ];
  state[4] = data[MHDTerm::BX];
  state[5] = data[MHDTerm::BY];
  state[6] = data[MHDTerm::BZ];
  const CFreal B = data[MHDTerm::B];
  const CFreal gammaMinus1 = getModel()->getGamma() - 1.0;
  state[7] = data[MHDTerm::P]/gammaMinus1 + 0.5*
    (rho*data[MHDTerm::V]*data[MHDTerm::V] + B*B);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD
