#include "Environment/ObjectProvider.hh"
#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DCons, ConvectiveVarSet, NavierStokesModule, 1>
euler2DConsProvider("Euler2DCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DCons::Euler2DCons(Common::SafePtr<BaseTerm> term) :
  Euler2DVarSet(term),
  _rightEv(4,4),
  _leftEv(4,4)
{
  vector<std::string> names(4);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoE";

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DCons::~Euler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setup()
{
  Euler2DVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setConstJacob()
{
  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  (*jacobians)[0](0,1) = 1.0;
  (*jacobians)[0](1,3) = gammaMinus1;

  (*jacobians)[1](0,2) = 1.0;
  (*jacobians)[1](2,3) = gammaMinus1;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal avH = linearData[EulerTerm::H];
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus3 = gamma - 3.;
  const CFreal uu = avU*avU;
  const CFreal uv = avU*avV;
  const CFreal vv = avV*avV;
  const CFreal sumVel2 = 0.5*gammaMinus1*(uu + vv);
  const CFreal vn = avU*nx + avV*ny;

  jacob(0,1) = nx;
  jacob(0,2) = ny;
  
  jacob(1,0) = (sumVel2 - uu)*nx - uv*ny;
  jacob(1,1) = -gammaMinus3*avU*nx + avV*ny;
  jacob(1,2) = -gammaMinus1*avV*nx + avU*ny;
  jacob(1,3) = gammaMinus1*nx;
  
  jacob(2,0) = -uv*nx + (sumVel2 - vv)*ny;
  jacob(2,1) = avV*nx - gammaMinus1*avU*ny;
  jacob(2,2) = avU*nx - gammaMinus3*avV*ny;
  jacob(2,3) = gammaMinus1*ny;
  
  jacob(3,0) = (sumVel2- avH)*vn;
  jacob(3,1) = (-gammaMinus1*uu + avH)*nx - gammaMinus1*uv*ny;
  jacob(3,2) = -gammaMinus1*uv*nx + (-gammaMinus1*vv + avH)*ny;
  jacob(3,3) = gamma*vn;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal avH = linearData[EulerTerm::H];
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus3 = gamma - 3.;
  const CFreal uu = avU*avU;
  const CFreal uv = avU*avV;
  const CFreal vv = avV*avV;
  const CFreal sumVel2 = 0.5*gammaMinus1*(uu + vv);

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](1,0) = sumVel2 - uu;
  (*jacobians)[0](1,1) = -gammaMinus3*avU;
  (*jacobians)[0](1,2) = -gammaMinus1*avV;
  (*jacobians)[0](2,0) = -uv;
  (*jacobians)[0](2,1) = avV;
  (*jacobians)[0](2,2) = avU;
  (*jacobians)[0](3,0) = sumVel2*avU - avU*avH;
  (*jacobians)[0](3,1) = -gammaMinus1*uu + avH;
  (*jacobians)[0](3,2) = -gammaMinus1*uv;
  (*jacobians)[0](3,3) = getModel()->getGamma()*avU;

  (*jacobians)[1](1,0) = -uv;
  (*jacobians)[1](1,1) = avV;
  (*jacobians)[1](1,2) = avU;
  (*jacobians)[1](2,0) = sumVel2 - vv;
  (*jacobians)[1](2,1) = -gammaMinus1*avU;
  (*jacobians)[1](2,2) = -gammaMinus3*avV;
  (*jacobians)[1](3,0) = sumVel2*avV - avV*avH;
  (*jacobians)[1](3,1) = -gammaMinus1*uv;
  (*jacobians)[1](3,2) = -gammaMinus1*vv + avH;
  (*jacobians)[1](3,3) = getModel()->getGamma()*avV;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avH   = linearData[EulerTerm::H];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal ovAvA = 1./avA;
  const CFreal ovAvA2 = ovAvA*ovAvA;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal um = avU*nx + avV*ny;
  const CFreal ra = 0.5*avRho*ovAvA;
  const CFreal coeffM2 = 0.5*gammaMinus1*(avU*avU + avV*avV)*ovAvA2;
  const CFreal ovAvRho = 1./avRho;
  const CFreal uDivA = gammaMinus1*avU*ovAvA;
  const CFreal vDivA = gammaMinus1*avV*ovAvA;
  const CFreal ovAvRhoA = ovAvRho*ovAvA;

  rightEv(0,0) = 1.;
  rightEv(0,1) = 0.;
  rightEv(0,2) = ra;
  rightEv(0,3) = ra;
  rightEv(1,0) = avU;
  rightEv(1,1) = avRho*ny;
  rightEv(1,2) = ra*(avU + avA*nx);
  rightEv(1,3) = ra*(avU - avA*nx);
  rightEv(2,0) = avV;
  rightEv(2,1) = -avRho*nx;;
  rightEv(2,2) = ra*(avV + avA*ny);
  rightEv(2,3) = ra*(avV - avA*ny);
  rightEv(3,0) = 0.5*(avU*avU +avV*avV);
  rightEv(3,1) = avRho*(avU*ny - avV*nx);
  rightEv(3,2) = ra*(avH + avA*um);
  rightEv(3,3) = ra*(avH - avA*um);

  leftEv(0,0) = 1.- coeffM2;
  leftEv(0,1) = uDivA*ovAvA;
  leftEv(0,2) = vDivA*ovAvA;
  leftEv(0,3) = -gammaMinus1*ovAvA2;
  leftEv(1,0) = ovAvRho*(avV*nx - avU*ny);
  leftEv(1,1) = ovAvRho*ny;
  leftEv(1,2) = -ovAvRho*nx;
  leftEv(1,3) = 0.0;
  leftEv(2,0) = avA*ovAvRho*(coeffM2 - um*ovAvA);
  leftEv(2,1) = ovAvRho*(nx - uDivA);
  leftEv(2,2) = ovAvRho*(ny - vDivA);
  leftEv(2,3) = gammaMinus1*ovAvRhoA;
  leftEv(3,0) = avA*ovAvRho*(coeffM2 + um*ovAvA);
  leftEv(3,1) = -ovAvRho*(nx + uDivA);
  leftEv(3,2) = -ovAvRho*(ny + vDivA);
  leftEv(3,3) = gammaMinus1*ovAvRhoA;
  
  eValues[0] = um;
  eValues[1] = um;
  eValues[2] = um + avA;
  eValues[3] = um - avA;

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
  << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
  << leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DChar::computeEigenValuesVectors" << "\n"
  << eValues << "\n" << "\n");
}


//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setEigenVect1(RealVector& r1,
				State& state,
				const RealVector& normal)
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];

  r1[0] = 1.0;
  r1[1] = u;
  r1[2] = v;
  r1[3] = 0.5*(u*u + v*v);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setEigenVect2(RealVector& r2,
				State& state,
				const RealVector& normal)
{
  r2[0] = 0.0;
  r2[1] = state[0]*normal[1];
  r2[2] = -state[0]*normal[0];
  r2[3] = state[1]*normal[1] - state[2]*normal[0];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setEigenVect3(RealVector& r3,
				State& state,
				const RealVector& normal)
{
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal rhoInv = 1./state[0];
  const CFreal u = state[1]*rhoInv;
  const CFreal v = state[2]*rhoInv;
  const CFreal uuvv =  0.5*(u*u + v*v);
  const CFreal p = gammaMinus1*(state[3] - state[0]*uuvv);
  const CFreal a = sqrt(gamma*p*rhoInv);
  const CFreal H = rhoInv*(gammaDivGammaMinus1*p + state[0]*uuvv);
  const CFreal ra = 0.5*state[0]/a;

  r3[0] = ra;
  r3[1] = ra*(u + a*normal[0]);
  r3[2] = ra*(v + a*normal[1]);
  r3[3] = ra*(H + a*(u*normal[0] + v*normal[1]));
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setEigenVect4(RealVector& r4,
				State& state,
				const RealVector& normal)
{
  const CFreal rhoInv = 1./state[0];
  const CFreal u = state[1]*rhoInv;
  const CFreal v = state[2]*rhoInv;
  const CFreal uuvv =  0.5*(u*u + v*v);

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;

  const CFreal p = gammaMinus1*(state[3] - state[0]*uuvv);
  const CFreal a = sqrt(gamma*p*rhoInv);
  const CFreal H = rhoInv*(gammaDivGammaMinus1*p + state[0]*uuvv);
  const CFreal ra = 0.5*state[0]/a;

  r4[0] = ra;
  r4[1] = ra*(u - a*normal[0]);
  r4[2] = ra*(v - a*normal[1]);
  r4[3] = ra*(H - a*(u*normal[0] + v*normal[1]));
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::splitJacobian(RealMatrix& jacobPlus,
				RealMatrix& jacobMin,
				RealVector& eValues,
				const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avH   = linearData[EulerTerm::H];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal um = avU*nx + avV*ny;
  const CFreal ra = 0.5*avRho/avA;
  const CFreal avA2 = avA*avA;
  const CFreal coeffM2 = 0.5*gammaMinus1 * (avU*avU + avV*avV)/avA2;
  const CFreal ovAvRho = 1./avRho;
  const CFreal uDivA = gammaMinus1*avU/avA;
  const CFreal vDivA = gammaMinus1*avV/avA;
  const CFreal rhoA = avRho*avA;

  _rightEv(0,0) = 1.;
  _rightEv(0,1) = 0.;
  _rightEv(0,2) = ra;
  _rightEv(0,3) = ra;
  _rightEv(1,0) = avU;
  _rightEv(1,1) = avRho*ny;
  _rightEv(1,2) = ra*(avU + avA*nx);
  _rightEv(1,3) = ra*(avU - avA*nx);
  _rightEv(2,0) = avV;
  _rightEv(2,1) = -avRho*nx;;
  _rightEv(2,2) = ra*(avV + avA*ny);
  _rightEv(2,3) = ra*(avV - avA*ny);
  _rightEv(3,0) = 0.5*(avU*avU +avV*avV);
  _rightEv(3,1) = avRho*(avU*ny - avV*nx);
  _rightEv(3,2) = ra*(avH + avA*um);
  _rightEv(3,3) = ra*(avH - avA*um);

  _leftEv(0,0) = 1.- coeffM2;
  _leftEv(0,1) = uDivA/avA;
  _leftEv(0,2) = vDivA/avA;
  _leftEv(0,3) = -gammaMinus1/avA2;
  _leftEv(1,0) = ovAvRho*(avV*nx - avU*ny);
  _leftEv(1,1) = ovAvRho*ny;
  _leftEv(1,2) = -ovAvRho*nx;
  _leftEv(1,3) = 0.0;
  _leftEv(2,0) = avA*ovAvRho*(coeffM2 - um/avA);
  _leftEv(2,1) = ovAvRho*(nx - uDivA);
  _leftEv(2,2) = ovAvRho*(ny - vDivA);
  _leftEv(2,3) = gammaMinus1/rhoA;
  _leftEv(3,0) = avA*ovAvRho*(coeffM2 + um/avA);
  _leftEv(3,1) = -ovAvRho*(nx + uDivA);
  _leftEv(3,2) = -ovAvRho*(ny + vDivA);
  _leftEv(3,3) = gammaMinus1/rhoA;

  eValues[0] = um;
  eValues[1] = um;
  eValues[2] = um + avA;
  eValues[3] = um - avA;
  // compute the eigen values + and -
  // for (CFuint iEq = 0; iEq < 4; ++iEq) {
  //     _eValuesP[iEq] = max(0.,eValues[iEq]);
  //     _eValuesM[iEq] = min(0.,eValues[iEq]);
  //   }

  // compute the eigen values + and -

// unused //  const CFreal speed = sqrt(avU*avU + avV*avV);
// unused //  const CFreal cosD = avU/speed;
// unused //  const CFreal sinD = avV/speed;
// unused //  const CFreal nEta = -nx*sinD + ny*cosD;

  // _jacobDissip = nEta;

  if (std::abs(_jacobDissip) > 0.0 || _delta > 0.0) {
    if (std::abs(_jacobDissip) > 0.0) {
      // modified eigenvalues to cure carbuncle
      const CFreal j2 = _jacobDissip*_jacobDissip;
      for (CFuint iEq = 0; iEq < eValues.size(); ++iEq) {
	_eValuesP[iEq] = max(0.,eValues[iEq]);
	const CFreal evP = _eValuesP[iEq];
	_eValuesP[iEq] = 0.5*(evP + sqrt(evP*evP + j2*avA2));

	_eValuesM[iEq] = min(0.,eValues[iEq]);
	const CFreal evM = _eValuesM[iEq];
	_eValuesM[iEq] = 0.5*(evM - sqrt(evM*evM + j2*avA2));
      }
    }
    else if (_delta > 0.0) {
      for (CFuint iEq = 0; iEq < eValues.size(); ++iEq) {
	const CFreal absEv = max(std::abs(eValues[iEq]), _delta);
	//  const CFreal absEv = std::abs(eValues[iEq]) + _delta;
	_eValuesP[iEq] = 0.5*(eValues[iEq] + absEv);
	_eValuesM[iEq] = 0.5*(eValues[iEq] - absEv);
      }
    }
  }
  else {
    for (CFuint iEq = 0; iEq < eValues.size(); ++iEq) {
      _eValuesP[iEq] = max(0.,eValues[iEq]);
      _eValuesM[iEq] = min(0.,eValues[iEq]);
    }
  }

  // compute jacobian + and -
  jacobPlus = _rightEv*(_eValuesP*_leftEv);
  jacobMin  = _rightEv*(_eValuesM*_leftEv);

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DCons::splitJacobian" << "\n"
  << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DCons::splitJacobian" << "\n"
  << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DCons::splitJacobian" << "\n"
  << eValues << "\n" << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::computePhysicalData(const State& state, RealVector& data)
{  
  // we assume that if conservative variables are used, the flow is compressible
  // p = static pressure
  const CFreal rho  = state[0];
  const CFreal ovRho = 1./rho;
  const CFreal u = state[1]*ovRho;
  const CFreal v = state[2]*ovRho;
  const CFreal V2 = u*u + v*v;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  
  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = gammaMinus1*(state[3] - 0.5*rho*V2);
  
  const CFreal pOvRho = data[EulerTerm::P]*ovRho;
  data[EulerTerm::E] = state[3]*ovRho;
  data[EulerTerm::H] = data[EulerTerm::E] + pOvRho;
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = pOvRho/getModel()->getR();
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::GAMMA] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
  const CFreal rho = data[EulerTerm::RHO];
  state[0] = rho;
  state[1] = rho*data[EulerTerm::VX];
  state[2] = rho*data[EulerTerm::VY];
  state[3] = rho*data[EulerTerm::H] - data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DCons::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  return sqrt(u*u + v*v);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::RHO];
  result[1] = state[1]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::RHO]*refData[EulerTerm::H];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/(refData[EulerTerm::RHO]);
  result[1] = state[1]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[2] = state[2]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[3] = state[3]/(refData[EulerTerm::RHO]*refData[EulerTerm::H]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCons::computePerturbedPhysicalData(const Framework::State& state,
					       const RealVector& pdataBkp,
					       RealVector& pdata,
					       CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool Euler2DCons::isValid(const RealVector& data)
{
  bool correct = true;
  enum index {RHO, RHOU,RHOV,RHOE};

  const CFreal rho = data[RHO];
  const CFreal ovRho = 1./rho;
  const CFreal u = data[RHOU]*ovRho;
  const CFreal v = data[RHOV]*ovRho;
  const CFreal V2 = u*u + v*v;

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  /// Next call returns R/M, in dimensional units.
  const CFreal R = getModel()->getR();

  const CFreal p = gammaMinus1*(data[RHOE] - 0.5*rho*V2);

  const CFreal T = p*ovRho/R;
  
  const CFreal a2 = gamma*p*ovRho;

  if( ( p < 0.) || ( T < 0.) || ( a2 < 0.) ){
  return correct = false;
  }

  return correct;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

