#include "MHD/MHD.hh"

#include "MHD3DProjectionCons.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshData.hh"
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

Environment::ObjectProvider<MHD3DProjectionCons, ConvectiveVarSet, MHDModule, 1> 
mhd3DProjectionConsProvider("MHD3DProjectionCons");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionCons::MHD3DProjectionCons(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DProjectionVarSet(term),
  _rm(),
  _rmInv(),
  _linearState(),
  _linState(),
  _stateFlux(),
  _dudw(),
  _rightMat()
{
  vector<std::string> names(9);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoW";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "rhoE";
  names[8] = "phi";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionCons::~MHD3DProjectionCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::setup()
{
  MHD3DProjectionVarSet::setup();

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

}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD3DProjectionCons::getExtraVarNames() const
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
  names[8] = "p";
  
  return names;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionCons::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                      RealMatrix& leftEv,
                                      RealVector& eValues,
                                      const RealVector& normal)
{
  // components of the rotated linearized state vector
  setRotationMatrices(normal, _rm, _rmInv);

  const RealVector& linearData = getModel()->getPhysicalData();

  //set the linear states in cons variables
  const CFreal rho = linearData[MHDProjectionTerm::RHO];
  _linearState[0] = rho;
  _linearState[1] = rho*linearData[MHDProjectionTerm::VX];
  _linearState[2] = rho*linearData[MHDProjectionTerm::VY];
  _linearState[3] = rho*linearData[MHDProjectionTerm::VZ];
  _linearState[4] = linearData[MHDProjectionTerm::BX];
  _linearState[5] = linearData[MHDProjectionTerm::BY];
  _linearState[6] = linearData[MHDProjectionTerm::BZ];
  const CFreal magB = linearData[MHDProjectionTerm::B];
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.);
  const CFreal rhoH = gammaDivGammaMinus1*linearData[MHDProjectionTerm::P]
    + 0.5*rho*linearData[MHDProjectionTerm::V]*linearData[MHDProjectionTerm::V] + 0.5*magB*magB;
  _linearState[7] = rhoH - linearData[MHDProjectionTerm::P];
  _linearState[8] = linearData[MHDProjectionTerm::PHI];

  _linState = _rm*_linearState;

  const CFreal gammaMinus1 = gamma - 1.;
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
  // const CFreal phiAv = _linState[8];

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

  const CFreal refSpeed = getModel()->getRefSpeed();

  // eigenvalues of the system

  eValues[0] = vnAv;
  eValues[1] = refSpeed;
  eValues[2] = -refSpeed;
  eValues[3] = vnAv + caAv;
  eValues[4] = vnAv - caAv;
  eValues[5] = vnAv + cfAv;
  eValues[6] = vnAv - cfAv;
  eValues[7] = vnAv + csAv;
  eValues[8] = vnAv - csAv;

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
  _rightMat(8,0) = 0.0;

  _rightMat(0,1) = 0.0;
  _rightMat(1,1) = 0.0;
  _rightMat(2,1) = 0.0;
  _rightMat(3,1) = 0.0;
  _rightMat(4,1) = 1.0;
  _rightMat(5,1) = 0.0;
  _rightMat(6,1) = 0.0;
  _rightMat(7,1) = 0.0;
  _rightMat(8,1) = refSpeed;

  _rightMat(0,2) = 0.0;
  _rightMat(1,2) = 0.0;
  _rightMat(2,2) = 0.0;
  _rightMat(3,2) = 0.0;
  _rightMat(4,2) = 1.0;
  _rightMat(5,2) = 0.0;
  _rightMat(6,2) = 0.0;
  _rightMat(7,2) = 0.0;
  _rightMat(8,2) = -refSpeed;

  _rightMat(0,3) = 0.0;
  _rightMat(1,3) = 0.0;
  _rightMat(2,3) = -bethaz*ovsq2;
  _rightMat(3,3) = bethap*ovsq2;
  _rightMat(4,3) = 0.0;
  _rightMat(5,3) = sqrt(rhoAv)*bethaz*ovsq2;
  _rightMat(6,3) = -sqrt(rhoAv)*bethap*ovsq2;
  _rightMat(7,3) = 0.0;
  _rightMat(8,3) = 0.0;

  _rightMat(0,4) = 0.0;
  _rightMat(1,4) = 0.0;
  _rightMat(2,4) = -bethaz*ovsq2;
  _rightMat(3,4) = bethap*ovsq2;
  _rightMat(4,4) = 0.0;
  _rightMat(5,4) = -sqrt(rhoAv)*bethaz*ovsq2;
  _rightMat(6,4) = sqrt(rhoAv)*bethap*ovsq2;
  _rightMat(7,4) = 0.0;
  _rightMat(8,4) = 0.0;

  _rightMat(0,5) = rhoAv*alphaf;
  _rightMat(1,5) = alphaf*cfAv;
  _rightMat(2,5) = -alphas*csAv*bethap*sgnbn;
  _rightMat(3,5) = -alphas*csAv*bethaz*sgnbn;
  _rightMat(4,5) = 0.0;
  _rightMat(5,5) = alphas*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,5) = alphas*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,5) = alphaf*gamma*pAv;
  _rightMat(8,5) = 0.0;

  _rightMat(0,6) = rhoAv*alphaf;
  _rightMat(1,6) = -alphaf*cfAv;
  _rightMat(2,6) = alphas*csAv*bethap*sgnbn;
  _rightMat(3,6) = alphas*csAv*bethaz*sgnbn;
  _rightMat(4,6) = 0.0;
  _rightMat(5,6) = alphas*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,6) = alphas*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,6) = alphaf*gamma*pAv;
  _rightMat(8,6) = 0.0;

  _rightMat(0,7) = rhoAv*alphas;
  _rightMat(1,7) = alphas*csAv;
  _rightMat(2,7) = alphaf*cfAv*bethap*sgnbn;
  _rightMat(3,7) = alphaf*cfAv*bethaz*sgnbn;
  _rightMat(4,7) = 0.0;
  _rightMat(5,7) = -alphaf*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,7) = -alphaf*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,7) = alphas*gamma*pAv;
  _rightMat(8,7) = 0.0;

  _rightMat(0,8) = rhoAv*alphas;
  _rightMat(1,8) = -alphas*csAv;
  _rightMat(2,8) = -alphaf*cfAv*bethap*sgnbn;
  _rightMat(3,8) = -alphaf*cfAv*bethaz*sgnbn;
  _rightMat(4,8) = 0.0;
  _rightMat(5,8) = -alphaf*sqrt(rhoAv)*aAv*bethap;
  _rightMat(6,8) = -alphaf*sqrt(rhoAv)*aAv*bethaz;
  _rightMat(7,8) = alphas*gamma*pAv;
  _rightMat(8,8) = 0.0;

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
  _dudw(8,0) = 0.0;

  _dudw(0,1) = 0.0;
  _dudw(1,1) = rhoAv;
  _dudw(2,1) = 0.0;
  _dudw(3,1) = 0.0;
  _dudw(4,1) = 0.0;
  _dudw(5,1) = 0.0;
  _dudw(6,1) = 0.0;
  _dudw(7,1) = rhoAv*vnAv;
  _dudw(8,1) = 0.0;

  _dudw(0,2) = 0.0;
  _dudw(1,2) = 0.0;
  _dudw(2,2) = rhoAv;
  _dudw(3,2) = 0.0;
  _dudw(4,2) = 0.0;
  _dudw(5,2) = 0.0;
  _dudw(6,2) = 0.0;
  _dudw(7,2) = rhoAv*vtAv;
  _dudw(8,2) = 0.0;

  _dudw(0,3) = 0.0;
  _dudw(1,3) = 0.0;
  _dudw(2,3) = 0.0;
  _dudw(3,3) = rhoAv;
  _dudw(4,3) = 0.0;
  _dudw(5,3) = 0.0;
  _dudw(6,3) = 0.0;
  _dudw(7,3) = rhoAv*vzAv;
  _dudw(8,3) = 0.0;

  _dudw(0,4) = 0.0;
  _dudw(1,4) = 0.0;
  _dudw(2,4) = 0.0;
  _dudw(3,4) = 0.0;
  _dudw(4,4) = 1.0;
  _dudw(5,4) = 0.0;
  _dudw(6,4) = 0.0;
  _dudw(7,4) = BnAv;
  _dudw(8,4) = 0.0;

  _dudw(0,5) = 0.0;
  _dudw(1,5) = 0.0;
  _dudw(2,5) = 0.0;
  _dudw(3,5) = 0.0;
  _dudw(4,5) = 0.0;
  _dudw(5,5) = 1.0;
  _dudw(6,5) = 0.0;
  _dudw(7,5) = BtAv;
  _dudw(8,5) = 0.0;

  _dudw(0,6) = 0.0;
  _dudw(1,6) = 0.0;
  _dudw(2,6) = 0.0;
  _dudw(3,6) = 0.0;
  _dudw(4,6) = 0.0;
  _dudw(5,6) = 0.0;
  _dudw(6,6) = 1.0;
  _dudw(7,6) = BzAv;
  _dudw(8,6) = 0.0;

  _dudw(0,7) = 0.0;
  _dudw(1,7) = 0.0;
  _dudw(2,7) = 0.0;
  _dudw(3,7) = 0.0;
  _dudw(4,7) = 0.0;
  _dudw(5,7) = 0.0;
  _dudw(6,7) = 0.0;
  _dudw(7,7) = 1.0/(gamma-1.0);
  _dudw(8,7) = 0.0;

  _dudw(0,8) = 0.0;
  _dudw(1,8) = 0.0;
  _dudw(2,8) = 0.0;
  _dudw(3,8) = 0.0;
  _dudw(4,8) = 0.0;
  _dudw(5,8) = 0.0;
  _dudw(6,8) = 0.0;
  _dudw(7,8) = 0.0;
  _dudw(8,8) = 1.0;

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
  leftEv(0,8) = 0.0;

  leftEv(1,0) = 0.0;
  leftEv(1,1) = 0.0;
  leftEv(1,2) = 0.0;
  leftEv(1,3) = 0.0;
  leftEv(1,4) = 1.0;
  leftEv(1,5) = 0.0;
  leftEv(1,6) = 0.0;
  leftEv(1,7) = 0.0;
  leftEv(1,8) = 1.0/refSpeed;

  leftEv(2,0) = 0.0;
  leftEv(2,1) = 0.0;
  leftEv(2,2) = 0.0;
  leftEv(2,3) = 0.0;
  leftEv(2,4) = 1.0;
  leftEv(2,5) = 0.0;
  leftEv(2,6) = 0.0;
  leftEv(2,7) = 0.0;
  leftEv(2,8) = -1.0/refSpeed;

  leftEv(3,0) = 0.0;
  leftEv(3,1) = 0.0;
  leftEv(3,2) = -bethaz*ovsq2;
  leftEv(3,3) = bethap*ovsq2;
  leftEv(3,4) = 0.0;
  leftEv(3,5) = bethaz*ovsq2/sqrt(rhoAv);
  leftEv(3,6) = -bethap*ovsq2/sqrt(rhoAv);
  leftEv(3,7) = 0.0;
  leftEv(3,8) = 0.0;

  leftEv(4,0) = 0.0;
  leftEv(4,1) = 0.0;
  leftEv(4,2) = -bethaz*ovsq2;
  leftEv(4,3) = bethap*ovsq2;
  leftEv(4,4) = 0.0;
  leftEv(4,5) = -bethaz*ovsq2/sqrt(rhoAv);
  leftEv(4,6) = bethap*ovsq2/sqrt(rhoAv);
  leftEv(4,7) = 0.0;
  leftEv(4,8) = 0.0;

  leftEv(5,0) = 0.0;
  leftEv(5,1) = alphaf*cfAv/(2.0*aAv2);
  leftEv(5,2) = -alphas*csAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(5,3) = -alphas*csAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(5,4) = 0.0;
  leftEv(5,5) = alphas*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(5,6) = alphas*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(5,7) = alphaf/(2.0*rhoAv*aAv2);
  leftEv(5,8) = 0.0;

  leftEv(6,0) = 0.0;
  leftEv(6,1) = -alphaf*cfAv/(2.0*aAv2);
  leftEv(6,2) = alphas*csAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(6,3) = alphas*csAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(6,4) = 0.0;
  leftEv(6,5) = alphas*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(6,6) = alphas*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(6,7) = alphaf/(2.0*rhoAv*aAv2);
  leftEv(6,8) = 0.0;

  leftEv(7,0) = 0.0;
  leftEv(7,1) = alphas*csAv/(2.0*aAv2);
  leftEv(7,2) = alphaf*cfAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(7,3) = alphaf*cfAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(7,4) = 0.0;
  leftEv(7,5) = -alphaf*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(7,6) = -alphaf*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(7,7) = alphas/(2.0*rhoAv*aAv2);
  leftEv(7,8) = 0.0;

  leftEv(8,0) = 0.0;
  leftEv(8,1) = -alphas*csAv/(2.0*aAv2);
  leftEv(8,2) = -alphaf*cfAv*bethap*sgnbn/(2.0*aAv2);
  leftEv(8,3) = -alphaf*cfAv*bethaz*sgnbn/(2.0*aAv2);
  leftEv(8,4) = 0.0;
  leftEv(8,5) = -alphaf*bethap/(2.0*sqrt(rhoAv)*aAv);
  leftEv(8,6) = -alphaf*bethaz/(2.0*sqrt(rhoAv)*aAv);
  leftEv(8,7) = alphas/(2.0*rhoAv*aAv2);
  leftEv(8,8) = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

CFuint MHD3DProjectionCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::computePhysicalData(const State& state,
					      RealVector& data)
{
  const CFreal rho = state[0];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal w = state[3]/rho;
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal phi = state[8];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;

  data[MHDProjectionTerm::VX] = u;
  data[MHDProjectionTerm::VY] = v;
  data[MHDProjectionTerm::VZ] = w;
  data[MHDProjectionTerm::BX] = Bx;
  data[MHDProjectionTerm::BY] = By;
  data[MHDProjectionTerm::BZ] = Bz;
  data[MHDProjectionTerm::B] = sqrt(B2);
  data[MHDProjectionTerm::RHO] = rho;
  data[MHDProjectionTerm::V] = sqrt(V2);
  data[MHDProjectionTerm::P] = (getModel()->getGamma() - 1.)*
    (state[7] - 0.5*(rho*V2 + B2));
  data[MHDProjectionTerm::PHI] = phi;
  data[MHDProjectionTerm::A] = sqrt(getModel()->getGamma()*data[MHDProjectionTerm::P]/rho);
  
  data[MHDProjectionTerm::GAMMA] = getModel()->getGamma();
  
  const RealVector& node = state.getCoordinates();
  data[MHDTerm::XP] =  node[XX];
  data[MHDTerm::YP] =  node[YY];
  data[MHDTerm::ZP] =  node[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::setRotationMatrices(const RealVector& normals,
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
  rm(8,8) = 1.;


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
  rmInv(8,8) = 1.;

  // Alternative rotation matrix and its inverse are commented out (this doesn't work 
  // with Roe scheme), it is kept for a comparison only in this function with 
  // the above rotation matrix used elsewhere in COOLFluiD and which works fine

// const CFreal nx = normals[0];
//   const CFreal ny = normals[1];
//   const CFreal nz = normals[2];

//   const CFreal cosThetaY = sqrt(nx*nx + ny*ny) / sqrt(nx*nx + ny*ny + nz*nz);
//   const CFreal sinThetaY = nz / sqrt(nx*nx + ny*ny + nz*nz);
//   const CFreal cosThetaZ = nx / sqrt(nx*nx + ny*ny);
//   const CFreal sinThetaZ = ny / sqrt(nx*nx + ny*ny);

//   rm(0,0) = 1.0;
//   rm(0,1) = 0.0;
//   rm(0,2) = 0.0;
//   rm(0,3) = 0.0;
//   rm(0,4) = 0.0;
//   rm(0,5) = 0.0;
//   rm(0,6) = 0.0;
//   rm(0,7) = 0.0;
//   rm(0,8) = 0.0;

//   rm(1,0) = 0.0;
//   rm(1,1) = cosThetaY*cosThetaZ;
//   rm(1,2) = cosThetaY*sinThetaZ;
//   rm(1,3) = sinThetaY;
//   rm(1,4) = 0.0;
//   rm(1,5) = 0.0;
//   rm(1,6) = 0.0;
//   rm(1,7) = 0.0;
//   rm(1,8) = 0.0;

//   rm(2,0) = 0.0;
//   rm(2,1) = -sinThetaZ;
//   rm(2,2) = cosThetaZ;
//   rm(2,3) = 0.0;
//   rm(2,4) = 0.0;
//   rm(2,5) = 0.0;
//   rm(2,6) = 0.0;
//   rm(2,7) = 0.0;
//   rm(2,8) = 0.0;

//   rm(3,0) = 0.0;
//   rm(3,1) = -sinThetaY*cosThetaZ;
//   rm(3,2) = -sinThetaY*sinThetaZ;
//   rm(3,3) = cosThetaY;
//   rm(3,4) = 0.0;
//   rm(3,5) = 0.0;
//   rm(3,6) = 0.0;
//   rm(3,7) = 0.0;
//   rm(3,8) = 0.0;

//   rm(4,0) = 0.0;
//   rm(4,1) = 0.0;
//   rm(4,2) = 0.0;
//   rm(4,3) = 0.0;
//   rm(4,4) = cosThetaY*cosThetaZ;
//   rm(4,5) = cosThetaY*sinThetaZ;
//   rm(4,6) = sinThetaY;
//   rm(4,7) = 0.0;
//   rm(4,8) = 0.0;

//   rm(5,0) = 0.0;
//   rm(5,1) = 0.0;
//   rm(5,2) = 0.0;
//   rm(5,3) = 0.0;
//   rm(5,4) = -sinThetaZ;
//   rm(5,5) = cosThetaZ;
//   rm(5,6) = 0.0;
//   rm(5,7) = 0.0;
//   rm(5,8) = 0.0;

//   rm(6,0) = 0.0;
//   rm(6,1) = 0.0;
//   rm(6,2) = 0.0;
//   rm(6,3) = 0.0;
//   rm(6,4) = -sinThetaY*cosThetaZ;
//   rm(6,5) = -sinThetaY*sinThetaZ;
//   rm(6,6) = cosThetaY;
//   rm(6,7) = 0.0;
//   rm(6,8) = 0.0;

//   rm(7,0) = 0.0;
//   rm(7,1) = 0.0;
//   rm(7,2) = 0.0;
//   rm(7,3) = 0.0;
//   rm(7,4) = 0.0;
//   rm(7,5) = 0.0;
//   rm(7,6) = 0.0;
//   rm(7,7) = 1.0;
//   rm(7,8) = 0.0;

//   rm(8,0) = 0.0;
//   rm(8,1) = 0.0;
//   rm(8,2) = 0.0;
//   rm(8,3) = 0.0;
//   rm(8,4) = 0.0;
//   rm(8,5) = 0.0;
//   rm(8,6) = 0.0;
//   rm(8,7) = 0.0;
//   rm(8,8) = 1.0;

//   rmInv(0,0) = 1.0;
//   rmInv(0,1) = 0.0;
//   rmInv(0,2) = 0.0;
//   rmInv(0,3) = 0.0;
//   rmInv(0,4) = 0.0;
//   rmInv(0,5) = 0.0;
//   rmInv(0,6) = 0.0;
//   rmInv(0,7) = 0.0;
//   rmInv(0,8) = 0.0;

//   rmInv(1,0) = 0.0;
//   rmInv(1,1) = cosThetaY*cosThetaZ;
//   rmInv(1,2) = -sinThetaZ;
//   rmInv(1,3) = -sinThetaY*cosThetaZ;
//   rmInv(1,4) = 0.0;
//   rmInv(1,5) = 0.0;
//   rmInv(1,6) = 0.0;
//   rmInv(1,7) = 0.0;
//   rmInv(1,8) = 0.0;

//   rmInv(2,0) = 0.0;
//   rmInv(2,1) = cosThetaY*sinThetaZ;
//   rmInv(2,2) = cosThetaZ;
//   rmInv(2,3) = -sinThetaY*sinThetaZ;
//   rmInv(2,4) = 0.0;
//   rmInv(2,5) = 0.0;
//   rmInv(2,6) = 0.0;
//   rmInv(2,7) = 0.0;
//   rmInv(2,8) = 0.0;

//   rmInv(3,0) = 0.0;
//   rmInv(3,1) = sinThetaY;
//   rmInv(3,2) = 0.0;
//   rmInv(3,3) = cosThetaY;
//   rmInv(3,4) = 0.0;
//   rmInv(3,5) = 0.0;
//   rmInv(3,6) = 0.0;
//   rmInv(3,7) = 0.0;
//   rmInv(3,8) = 0.0;

//   rmInv(4,0) = 0.0;
//   rmInv(4,1) = 0.0;
//   rmInv(4,2) = 0.0;
//   rmInv(4,3) = 0.0;
//   rmInv(4,4) = cosThetaY*cosThetaZ;
//   rmInv(4,5) = -sinThetaZ;
//   rmInv(4,6) = -sinThetaY*cosThetaZ;
//   rmInv(4,7) = 0.0;
//   rmInv(4,8) = 0.0;

//   rmInv(5,0) = 0.0;
//   rmInv(5,1) = 0.0;
//   rmInv(5,2) = 0.0;
//   rmInv(5,3) = 0.0;
//   rmInv(5,4) = cosThetaY*sinThetaZ;
//   rmInv(5,5) = cosThetaZ;
//   rmInv(5,6) = -sinThetaY*sinThetaZ;
//   rmInv(5,7) = 0.0;
//   rmInv(5,8) = 0.0;

//   rmInv(6,0) = 0.0;
//   rmInv(6,1) = 0.0;
//   rmInv(6,2) = 0.0;
//   rmInv(6,3) = 0.0;
//   rmInv(6,4) = sinThetaY;
//   rmInv(6,5) = 0.0;
//   rmInv(6,6) = cosThetaY;
//   rmInv(6,7) = 0.0;
//   rmInv(6,8) = 0.0;

//   rmInv(7,0) = 0.0;
//   rmInv(7,1) = 0.0;
//   rmInv(7,2) = 0.0;
//   rmInv(7,3) = 0.0;
//   rmInv(7,4) = 0.0;
//   rmInv(7,5) = 0.0;
//   rmInv(7,6) = 0.0;
//   rmInv(7,7) = 1.0;
//   rmInv(7,8) = 0.0;

//   rmInv(8,0) = 0.0;
//   rmInv(8,1) = 0.0;
//   rmInv(8,2) = 0.0;
//   rmInv(8,3) = 0.0;
//   rmInv(8,4) = 0.0;
//   rmInv(8,5) = 0.0;
//   rmInv(8,6) = 0.0;
//   rmInv(8,7) = 0.0;
//   rmInv(8,8) = 1.0;

// // check for rm * rmInv = I

// //   RealMatrix matrixMultip(PhysicalModel::getNbEq(),
// //                        PhysicalModel::getNbEq());

// //   const CFuint nbEqs = PhysicalModel::getNbEq();

// //   matrixMultip = rm*rmInv;

// //   for(int i = 0; i < nbEqs; i++)
// //     for(int j = 0; j < nbEqs; j++)
// //       if (matrixMultip(i,j) != 1.0)
// //      cout << "here " << matrixMultip(i,j) << endl;

  CFLogDebugMax( "RM*RMI = " << "\n" << RealMatrix(rm*rmInv) << "\n");
  CFLogDebugMax( "RMI*RM = " << "\n" << RealMatrix(rmInv*rm) << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::setWaveStrengths(RealVector& waveStrengths,
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

  RealVector leftStateVector(Framework::PhysicalModelStack::getActive()->getNbEq());
  RealVector rightStateVector(Framework::PhysicalModelStack::getActive()->getNbEq());
  RealVector leftStateVectorRot(Framework::PhysicalModelStack::getActive()->getNbEq());
  RealVector rightStateVectorRot(Framework::PhysicalModelStack::getActive()->getNbEq());

  leftStateVector[0] = leftState[0];
  leftStateVector[1] = leftState[1];
  leftStateVector[2] = leftState[2];
  leftStateVector[3] = leftState[3];
  leftStateVector[4] = leftState[4];
  leftStateVector[5] = leftState[5];
  leftStateVector[6] = leftState[6];
  leftStateVector[7] = leftState[7];
  leftStateVector[8] = leftState[8];

  rightStateVector[0] = rightState[0];
  rightStateVector[1] = rightState[1];
  rightStateVector[2] = rightState[2];
  rightStateVector[3] = rightState[3];
  rightStateVector[4] = rightState[4];
  rightStateVector[5] = rightState[5];
  rightStateVector[6] = rightState[6];
  rightStateVector[7] = rightState[7];
  rightStateVector[8] = rightState[8];

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
  leftState[8] = leftStateVectorRot[8];

  rightState[0] = rightStateVectorRot[0];
  rightState[1] = rightStateVectorRot[1];
  rightState[2] = rightStateVectorRot[2];
  rightState[3] = rightStateVectorRot[3];
  rightState[4] = rightStateVectorRot[4];
  rightState[5] = rightStateVectorRot[5];
  rightState[6] = rightStateVectorRot[6];
  rightState[7] = rightStateVectorRot[7];
  rightState[8] = rightStateVectorRot[8];

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
  const CFreal phiL = leftState[8];
  const CFreal phiR = rightState[8];
  const CFreal v2L = vnL*vnL + vtL*vtL + vzL*vzL;
  const CFreal v2R = vnR*vnR + vtR*vtR + vzR*vzR;
  const CFreal B2L = BnL*BnL + BtL*BtL + BzL*BzL;
  const CFreal B2R = BnR*BnR + BtR*BtR + BzR*BzR;

  const CFreal pL = gammaMinus1*(rhoEL - 0.5*rhoL*v2L - 0.5*B2L);
  const CFreal pR = gammaMinus1*(rhoER - 0.5*rhoR*v2R - 0.5*B2R);

  RealVector dw(PhysicalModelStack::getActive()->getNbEq());

  dw[0] = rhoR - rhoL;
  dw[1] = vnR - vnL;
  dw[2] = vtR - vtL;
  dw[3] = vzR - vzL;
  dw[4] = BnR - BnL;
  dw[5] = BtR - BtL;
  dw[6] = BzR - BzL;
  dw[7] = pR - pL;
  dw[8] = phiR - phiL;

  waveStrengths = leftEv*dw;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::setDimensionalValuesPlusExtraValues(const State& state, 
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

  const CFreal B1dotB0 = state[4]*Bx0 + state[5]*By0 + state[6]*Bz0;
  const CFreal sqB0 = Bx0*Bx0 + By0*By0 + Bz0*Bz0;
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
  result[8] = state[8];

  extra.resize(10);
  
  extra[0] = Bx0;
  extra[1] = By0;
  extra[2] = Bz0;
  extra[3] = state[4] + extra[0];
  extra[4] = state[5] + extra[1];
  extra[5] = state[6] + extra[2];
  extra[6] = sqrt(extra[3]*extra[3] + extra[4]*extra[4] + extra[5]*extra[5]);
  extra[7] = state[7] + B1dotB0 + 0.5*sqB0;
					 
  extra[8] = (getModel()->getGamma() - 1.)*
    (state[7] - 0.5*(rho*V2 + sqB1));

  // std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  // std::string datahandleName = nsp + "_divBNodal";
  
  // DataHandle<CFreal> divBNodal = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);
  
  // const CFuint stateID = state.getLocalID();
  
  // extra[9] = divBNodal[stateID];

}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHD3DProjectionCons::getRotatedFlux(const Framework::State& state,
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
  const CFreal phiL = state[8];
  const CFreal v2L = vnL*vnL + vtL*vtL + vzL*vzL;
  const CFreal B2L = BnL*BnL + BtL*BtL + BzL*BzL;
  const CFreal pL = gammaMinus1*(rhoEL - 0.5*rhoL*v2L - 0.5*B2L);

  const CFreal refSpeed = getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  _stateFlux[0] = rhoL*vnL;
  _stateFlux[1] = rhoL*vnL*vnL + (pL+0.5*B2L) - BnL*BnL;
  _stateFlux[2] = rhoL*vnL*vtL - BnL*BtL;
  _stateFlux[3] = rhoL*vnL*vzL - BnL*BzL;
  _stateFlux[4] = phiL;
  _stateFlux[5] = vnL*BtL - vtL*BnL;
  _stateFlux[6] = vnL*BzL - vzL*BnL;
  _stateFlux[7] = (rhoEL+pL+0.5*B2L)*vnL-BnL*(vnL*BnL+vtL*BtL+vzL*BzL);
  _stateFlux[8] = refSpeedSq*BnL;

  return _stateFlux;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionCons::computeStateFromPhysicalData(const RealVector& data,
						  State& state)
{
  const CFreal rho = data[MHDProjectionTerm::RHO];
  const CFreal phi = data[MHDProjectionTerm::PHI];
  state[0] = rho;
  state[1] = rho*data[MHDProjectionTerm::VX];
  state[2] = rho*data[MHDProjectionTerm::VY];
  state[3] = rho*data[MHDProjectionTerm::VZ];
  state[4] = data[MHDProjectionTerm::BX];
  state[5] = data[MHDProjectionTerm::BY];
  state[6] = data[MHDProjectionTerm::BZ];
  const CFreal B = data[MHDProjectionTerm::B];
  state[7] = data[MHDProjectionTerm::P]/(getModel()->getGamma() - 1.0) + 0.5*
    (rho*data[MHDProjectionTerm::V]*data[MHDProjectionTerm::V] + B*B);
  state[8] = phi;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD
