#include "Environment/ObjectProvider.hh"
#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/Euler1DCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler1DCons, ConvectiveVarSet, NavierStokesModule, 1>
euler1DConsProvider("Euler1DCons");

//////////////////////////////////////////////////////////////////////////////

Euler1DCons::Euler1DCons(Common::SafePtr<BaseTerm> term) :
  Euler1DVarSet(term),
  _rightEv(3,3),
  _leftEv(3,3)
{
  vector<std::string> names(3);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoE";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler1DCons::~Euler1DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::setup()
{
  Euler1DVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::setConstJacob()
{
   vector<RealMatrix>* const jacobians =
      PhysicalModelStack::getActive()->getImplementor()->getJacobians();
  
   const CFreal gamma = getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.;
  
   (*jacobians)[0](0,1) = 1.0;
   (*jacobians)[0](1,2) = gammaMinus1;

}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
   const RealVector& linearData = getModel()->getPhysicalData();
   const CFreal nx = normal[XX];
   const CFreal avU = linearData[EulerTerm::VX];
   const CFreal avH = linearData[EulerTerm::H];
   const CFreal gamma = getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.;
   const CFreal gammaMinus3 = gamma - 3.;
   const CFreal uu = avU*avU;
   const CFreal sumVel2 = 0.5*gammaMinus1*uu;
   const CFreal vn = avU*nx;

   jacob(0,0) = 0.;
   jacob(0,1) = nx;
   jacob(0,2) = 0.;
   jacob(1,0) = 0.5*gammaMinus3*uu*nx;
   jacob(1,1) = -gammaMinus3*vn;
   jacob(1,2) = gammaMinus1*nx;
   jacob(2,0) = (sumVel2 - avH)*vn;
   jacob(2,1) = (avH - gammaMinus1*uu)*nx;
   jacob(2,2) = gamma*vn;
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::computeJacobians()
{
   const RealVector& linearData = getModel()->getPhysicalData();

   const CFreal avU = linearData[EulerTerm::VX];
   const CFreal avH = linearData[EulerTerm::H];
   const CFreal gamma = getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.;
   const CFreal gammaMinus3 = gamma - 3.;
   const CFreal uu = avU*avU;
   const CFreal sumVel2 = 0.5*gammaMinus1*uu;

   vector<RealMatrix>* const jacobians =
     PhysicalModelStack::getActive()->getImplementor()->getJacobians();

   (*jacobians)[0](1,0) = 0.5*gammaMinus3*uu;
   (*jacobians)[0](1,1) = -gammaMinus3*avU;
   (*jacobians)[0](2,0) = avU*(sumVel2 - avH);
   (*jacobians)[0](2,1) = avH - gammaMinus1*uu;
   (*jacobians)[0](2,2) = gamma*avU;

}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
   const RealVector& linearData = getModel()->getPhysicalData();

   const CFreal nx = normal[XX];
   const CFreal avRho = linearData[EulerTerm::RHO];
   const CFreal avU   = linearData[EulerTerm::VX];
   const CFreal avH   = linearData[EulerTerm::H];
   const CFreal avA   = linearData[EulerTerm::A];

   const CFreal gamma = getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.;
   
   const CFreal um = avU*nx;
   const CFreal ra = 0.5*avRho/avA;
   const CFreal avA2 = avA*avA;
   const CFreal coeffM2 = 0.5*gammaMinus1*avU*avU/avA2;
   const CFreal invAvRho = 1./avRho;
   const CFreal uDivA = gammaMinus1*avU/avA;
   const CFreal rhoA = avRho*avA;

   rightEv(0,0) = 1.;
   rightEv(0,1) = ra;
   rightEv(0,2) = ra;
   rightEv(1,0) = avU;
   rightEv(1,1) = ra*(avU + avA*nx);
   rightEv(1,2) = ra*(avU - avA*nx);
   rightEv(2,0) = 0.5*avU*avU;
   rightEv(2,1) = ra*(avH + avA*um);
   rightEv(2,2) = ra*(avH - avA*um);

   leftEv(0,0) = 1.- coeffM2;
   leftEv(0,1) = uDivA/avA;
   leftEv(0,2) = -gammaMinus1/avA2;
   leftEv(1,0) = avA*invAvRho*(coeffM2 - um/avA);
   leftEv(1,1) = invAvRho*(nx - uDivA);
   leftEv(1,2) = gammaMinus1/rhoA;
   leftEv(2,0) = avA*invAvRho*(coeffM2 + um/avA);
   leftEv(2,1) = -invAvRho*(nx + uDivA);
   leftEv(2,2) = gammaMinus1/rhoA;

   eValues[0] = um;
   eValues[1] = um + avA;
   eValues[2] = um - avA;
  
   //degugging infos
   CFLogDebugMax( "RightEigenvectors @Euler1DChar::computeEigenValuesVectors" << "\n"
   << rightEv << "\n");
   CFLogDebugMax( "LeftEigenvectors @Euler1DChar::computeEigenValuesVectors" << "\n"
   << leftEv << "\n");
   CFLogDebugMax( "EigenValues @Euler1DChar::computeEigenValuesVectors" << "\n"
   << eValues << "\n" << "\n");
}


//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::setEigenVect1(RealVector& r1,
        State& state,
        const RealVector& normal)
{
  const CFreal u = state[1]/state[0];

  r1[0] = 1.0;
  r1[1] = u;
  r1[2] = 0.5*u*u;
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::setEigenVect2(RealVector& r2,
        State& state,
        const RealVector& normal)
{
  const CFreal rhoInv = 1./state[0];
  const CFreal u = state[1]*rhoInv;
  const CFreal uu =  0.5*u*u;

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;

  const CFreal p = gammaMinus1*(state[2] - state[0]*uu);
  const CFreal a = sqrt(gamma*p*rhoInv);
  const CFreal H = rhoInv*(gammaDivGammaMinus1*p + state[0]*uu);
  const CFreal ra = 0.5*state[0]/a;

  r2[0] = ra;
  r2[1] = ra*(u + a*normal[0]);
  r2[2] = ra*(H + a*u*normal[0]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::setEigenVect3(RealVector& r3,
        State& state,
        const RealVector& normal)
{
  const CFreal rhoInv = 1./state[0];
  const CFreal u = state[1]*rhoInv;
  const CFreal uu =  0.5*u*u;

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;

  const CFreal p = gammaMinus1*(state[2] - state[0]*uu);
  const CFreal a = sqrt(gamma*p*rhoInv);
  const CFreal H = rhoInv*(gammaDivGammaMinus1*p + state[0]*uu);
  const CFreal ra = 0.5*state[0]/a;

  r3[0] = ra;
  r3[1] = ra*(u - a*normal[0]);
  r3[2] = ra*(H - a*u*normal[0]);
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler1DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::splitJacobian(RealMatrix& jacobPlus,
				RealMatrix& jacobMin,
				RealVector& eValues,
				const RealVector& normal)
{
   const RealVector& linearData = getModel()->getPhysicalData();

   const CFreal gamma = getModel()->getGamma();
   const CFreal gammaMinus1 = gamma - 1.;
   const CFreal nx = normal[XX];
   const CFreal avRho = linearData[EulerTerm::RHO];
   const CFreal avU   = linearData[EulerTerm::VX];
   const CFreal avH   = linearData[EulerTerm::H];
   const CFreal avA   = linearData[EulerTerm::A];
   const CFreal um = avU*nx;
   const CFreal ra = 0.5*avRho/avA;
   const CFreal avA2 = avA*avA;
   const CFreal coeffM2 = 0.5*gammaMinus1*avU*avU/avA2;
   const CFreal invAvRho = 1./avRho;
   const CFreal uDivA = gammaMinus1*avU/avA;
   const CFreal rhoA = avRho*avA;

   _rightEv(0,0) = 1.;
   _rightEv(0,1) = ra;
   _rightEv(0,2) = ra;
   _rightEv(1,0) = avU;
   _rightEv(1,1) = ra*(avU + avA*nx);
   _rightEv(1,2) = ra*(avU - avA*nx);
   _rightEv(2,0) = 0.5*avU*avU;
   _rightEv(2,1) = ra*(avH + avA*um);
   _rightEv(2,2) = ra*(avH - avA*um);

   _leftEv(0,0) = 1.- coeffM2;
   _leftEv(0,1) = uDivA/avA;
   _leftEv(0,2) = -gammaMinus1/avA2;
   _leftEv(1,0) = avA*invAvRho*(coeffM2 - um/avA);
   _leftEv(1,1) = invAvRho*(nx - uDivA);
   _leftEv(1,2) = gammaMinus1/rhoA;
   _leftEv(2,0) = avA*invAvRho*(coeffM2 + um/avA);
   _leftEv(2,1) = -invAvRho*(nx + uDivA);
   _leftEv(2,2) = gammaMinus1/rhoA;

   eValues[0] = um;
   eValues[1] = um + avA;
   eValues[2] = um - avA;
//   // compute the eigen values + and -
//   // for (CFuint iEq = 0; iEq < 4; ++iEq) {
//   //     _eValuesP[iEq] = max(0.,eValues[iEq]);
//   //     _eValuesM[iEq] = min(0.,eValues[iEq]);
//   //   }

//   // compute the eigen values + and -

// // unused //  const CFreal speed = sqrt(avU*avU + avV*avV);
// // unused //  const CFreal cosD = avU/speed;
// // unused //  const CFreal sinD = avV/speed;
// // unused //  const CFreal nEta = -nx*sinD + ny*cosD;

//   // _jacobDissip = nEta;

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

//   // compute jacobian + and -
   jacobPlus = _rightEv*(_eValuesP*_leftEv);
   jacobMin  = _rightEv*(_eValuesM*_leftEv);

   //degugging infos
   CFLogDebugMax( "RightEigenvectors @Euler1DCons::splitJacobian" << "\n"
   << _rightEv << "\n");
   CFLogDebugMax( "LeftEigenvectors @Euler1DCons::splitJacobian" << "\n"
   << _leftEv << "\n");
   CFLogDebugMax( "EigenValues @Euler1DCons::splitJacobian" << "\n"
   << eValues << "\n" << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::computePhysicalData(const State& state, RealVector& data)
{
  const CFreal rho  = state[0];
  const CFreal ovRho = 1./rho;
  const CFreal u = state[1]*ovRho;
  const CFreal V2 = u*u;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
    
  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = gammaMinus1*(state[2] - 0.5*rho*V2);
   
  const CFreal pOvRho = data[EulerTerm::P]*ovRho;
  
  data[EulerTerm::E] = state[2]*ovRho;
  data[EulerTerm::H] = data[EulerTerm::E] + pOvRho;
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = pOvRho/getModel()->getR();
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::GAMMA] = gamma;
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::computeStateFromPhysicalData(const RealVector& data, State& state)
{
  const CFreal rho = data[EulerTerm::RHO];
  state[0] = rho;
  state[1] = rho*data[EulerTerm::VX];
  state[2] = rho*data[EulerTerm::H] - data[EulerTerm::P];
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFreal Euler1DCons::getSpeed(const State& state) const
{
  return state[1]/state[0];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  
  result[0] = state[0]*refData[EulerTerm::RHO];
  result[1] = state[1]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::RHO]*refData[EulerTerm::H];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/(refData[EulerTerm::RHO]);
  result[1] = state[1]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[2] = state[2]/(refData[EulerTerm::RHO]*refData[EulerTerm::H]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
