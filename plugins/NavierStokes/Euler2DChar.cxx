#include "NavierStokes/NavierStokes.hh"
#include <numeric>

#include "Euler2DChar.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathConsts.hh"
#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DChar, ConvectiveVarSet, NavierStokesModule, 1>
euler2DCharProvider("Euler2DChar");

//////////////////////////////////////////////////////////////////////////////

Euler2DChar::Euler2DChar(Common::SafePtr<BaseTerm> term) :
  Euler2DVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DChar::~Euler2DChar()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DChar::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal speed = sqrt(avU*avU + avV+avV);
  const CFreal M     = speed/avA;
  const CFreal M2    = M*M;
  const CFreal eps    = 0.05; // this should be an argument of the constructor
  const CFreal beta   = sqrt(max(eps*eps,std::abs(M2 - 1.)));
  const CFreal beta2  = beta*beta;
  const CFreal chi    = beta/max(M, 1.);
  const CFreal nuPlus = 0.5*(M2 - 1. + beta2)/beta2;
  CFreal nuMin  = 0.5*(M2 - 1. - beta2)/beta2;

  vector<RealMatrix>* const jacobians = PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](0,0) = speed*chi*nuPlus;
  (*jacobians)[0](0,1) = speed*chi*nuMin;
  (*jacobians)[0](1,0) = speed*chi*nuMin;
  (*jacobians)[0](1,1) = speed*chi*nuPlus;
  (*jacobians)[0](2,2) = speed;
  (*jacobians)[0](3,3) = speed;

  (*jacobians)[1](0,0) = speed*chi/beta;
  (*jacobians)[1](1,1) = -speed*chi/beta;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DChar::computeEigenValuesVectors(RealMatrix& rightEv,
					RealMatrix& leftEv,
					RealVector& eValues,
					const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU  = linearData[EulerTerm::VX];
  const CFreal avV  = linearData[EulerTerm::VY];
  const CFreal avA  = linearData[EulerTerm::A];
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal speed = sqrt(avU*avU + avV*avV);
  const CFreal cosD = avU/speed;
  const CFreal sinD = avV/speed;
  const CFreal nCsi = nx*cosD + ny*sinD;
  const CFreal nEta = -nx*sinD + ny*cosD;
  const CFreal M     = speed/avA;
  const CFreal M2    = M*M;
  const CFreal eps   = 0.05; // this should be an argument of the constructor
  const CFreal beta   = sqrt(max<const CFreal>(eps*eps,std::abs(M2 - 1.)));
  const CFreal beta2  = beta*beta;
  const CFreal chi    = beta/max<const CFreal>(M,1.);
  CFreal nuPlus = 0.5*(M2 - 1. + beta2)/beta2;
  CFreal nuMin  = 0.5*(M2 - 1. - beta2)/beta2;
  if (M < 1.0) {
    nuMin = -1.0;
    nuPlus = 0.0;
  }
  const CFreal coeff  = sqrt(nuMin*nuMin*nCsi*nCsi + nEta*nEta/beta2);
  const CFreal um = speed*nCsi;
  const CFuint nbCoupledEq = 2;

  if (M > 1.0) { //if flow is fully supersonic, matrix is fully diagonalized
    for (CFuint i = 0; i < nbCoupledEq; ++i) {
      for (CFuint j = 0; j < nbCoupledEq; ++j) {
	rightEv(i,j) = (i == j) ? 1.0 : 0.0;
	leftEv(i,j) = (i == j) ? 1.0 : 0.0;
      }
    }

    // eValues include case with M = 1

    eValues[0] = speed*chi*(nCsi + nEta/beta);
    eValues[1] = speed*chi*(nCsi - nEta/beta);
  }
  else {//if flow is subsonic, matrix is fully diagonalized only if nCsi = 0
    if (std::abs(nCsi) > MathTools::MathConsts::CFrealEps()) {
      const CFreal k1 = -(nuMin*nCsi)/(nEta/beta - coeff);
      const CFreal k2 = -(nEta/beta + coeff)/(nuMin*nCsi);
      const CFreal det2 = k1*k2 - 1.0;

      rightEv(0,0) = k1;
      rightEv(0,1) = 1.0;
      rightEv(1,0) = 1.0;
      rightEv(1,1) = k2;

      leftEv(0,0) = k2/det2;
      leftEv(0,1) = -1.0/det2;
      leftEv(1,0) = -1.0/det2;
      leftEv(1,1) = k1/det2;
    }
    else {
      for (CFuint i = 0; i < nbCoupledEq; i++) {
	for (CFuint j = 0; j < nbCoupledEq; j++) {
	  rightEv(i,j) = (i == j) ? 1.0 : 0.0;
	  leftEv(i,j) = (i == j) ? 1.0 : 0.0;
	}
      }
    }
    eValues[0] = speed*chi*(nuPlus*nCsi + coeff);
    eValues[1] = speed*chi*(nuPlus*nCsi - coeff);
  }
  
  if (rightEv.nbRows() == 4) {
    rightEv(2,2) = rightEv(3,3) = 1.0;
    leftEv(2,2) = leftEv(3,3) = 1.0;
    eValues[2] = eValues[3] = um;
  }

#if 0  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
		 << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
	<< leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DChar::computeEigenValuesVectors" << "\n"
	<< eValues << "\n" << "\n");

  CFLog(DEBUGMAX, "The Matrix: " << "\n" <<
  	       leftEv*((*jacobians)[0]*nCsi +
  		       (*jacobians)[1]*nEta)*rightEv << "\n");

  CFLog(DEBUGMAX, "\n" << "Identity: " << "\n"
  	       << leftEv*rightEv << "\n"
  	       << "Identity: " << "\n" <<
  	       rightEv*leftEv << "\n");
#endif
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DChar::getBlockSeparator() const
{
  return 2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DChar::splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal)
{
  cf_assert(eValues.size() == 4);
  const RealVector& linearData = getModel()->getPhysicalData();
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal avU  = linearData[EulerTerm::VX];
  const CFreal avV  = linearData[EulerTerm::VY];
  const CFreal speed = linearData[EulerTerm::V];
  const CFreal cosD = avU/speed;
  const CFreal sinD = avV/speed;
  const CFreal nCsi = nx*cosD + ny*sinD;
  const CFreal nEta = -nx*sinD + ny*cosD;
  const CFreal M  = speed/linearData[EulerTerm::A];
  const CFreal M2 = M*M;
  const CFreal eps = 0.05;
  const CFreal beta   = sqrt(max(eps*eps,std::abs(M2 - 1.)));
  const CFreal beta2  = beta*beta;
  const CFreal chi    = beta/max(M,1.);
  CFreal nuPlus = 0.5*(M2 - 1. + beta2)/beta2;
  CFreal nuMin  = 0.5*(M2 - 1. - beta2)/beta2;
  if (M < 1.0) {
    nuMin = -1.0;
    nuPlus = 0.0;
  }
  const CFreal coeff  = sqrt(nuMin*nuMin*nCsi*nCsi + nEta*nEta/beta2);
  const CFreal um = speed*nCsi;

  if (MathTools::MathChecks::isZero(M-1.)) {
    CFout << "M = 1" << "\n";
  }

  // compute all the eigen values
  if (M > 1.0) {
    eValues[0] = speed*chi*(nCsi + nEta/beta);
    eValues[1] = speed*chi*(nCsi - nEta/beta);
  }
  else {
    eValues[0] = speed*chi*(nuPlus*nCsi + coeff);
    eValues[1] = speed*chi*(nuPlus*nCsi - coeff);
  }
  eValues[2] = um;
  eValues[3] = um;

  // compute the eigen values + and -
  for (CFuint iEq = 0; iEq < 4; ++iEq) {
    _eValuesP[iEq] = max(0.,eValues[iEq]);
    _eValuesM[iEq] = min(0.,eValues[iEq]);
  }

  // compute the "upwind" jacobians
  if (M > 1.0) {
    jacobPlus(0,0) = _eValuesP[0];
    jacobPlus(0,1) = 0.;
    jacobPlus(1,0) = 0.;
    jacobPlus(1,1) = _eValuesP[1];

    jacobMin(0,0) = _eValuesM[0];
    jacobMin(0,1) = 0.;
    jacobMin(1,0) = 0.;
    jacobMin(1,1) = _eValuesM[1];
  }
  else {
    // M <= 1.0 is M=1 considered here for the eigenvectors???
    if (std::abs(nCsi) > MathTools::MathConsts::CFrealEps()) {
      const CFreal k1 = -(nuMin*nCsi)/(nEta/beta - coeff);
      const CFreal k2 = -(nEta/beta + coeff)/(nuMin*nCsi);
      const CFreal k1k2 = k1*k2;
      const CFreal invDet2 = 1./(k1*k2 - 1.0);
      CFreal e1 = _eValuesP[0];
      CFreal e2 = _eValuesP[1];

      jacobPlus(0,0) = invDet2*(e1*k1k2 - e2);
      jacobPlus(0,1) = invDet2*k1*(-e1 + e2);
      jacobPlus(1,0) = invDet2*k2*(e1 - e2);
      jacobPlus(1,1) = invDet2*(-e1 + e2*k1k2);

      e1 = _eValuesM[0];
      e2 = _eValuesM[1];

      jacobMin(0,0) = invDet2*(e1*k1k2 - e2);
      jacobMin(0,1) = invDet2*k1*(-e1 + e2);
      jacobMin(1,0) = invDet2*k2*(e1 - e2);
      jacobMin(1,1) = invDet2*(-e1 + e2*k1k2);
    }
    else {
      // nCsi = 0 => nx*u + ny*v = 0 => Vn = 0
      CFLogDebugMax( "case nCsi = 0" << "\n");

      jacobPlus(0,0) = _eValuesP[0];
      jacobPlus(0,1) = 0.;
      jacobPlus(1,0) = 0.;
      jacobPlus(1,1) = _eValuesP[1];

      jacobMin(0,0) = _eValuesM[0];
      jacobMin(0,1) = 0.;
      jacobMin(1,0) = 0.;
      jacobMin(1,1) = _eValuesM[1];
    }
  }

  if (jacobPlus.nbRows() == 4) {
    jacobPlus(2,2) = _eValuesP[2];
    jacobPlus(3,3) = _eValuesP[3];

    jacobMin(2,2) = _eValuesM[2];
    jacobMin(3,3) = _eValuesM[3];
  }

#if 0 // this is added for testing purposes

  computeJacobians();

  RealMatrix m1(4,4);
  RealMatrix m2(4,4);
  RealMatrix m3(4,4);
  RealMatrix rightEv(4,4);
  RealMatrix leftEv(4,4);

  vector<RealMatrix>* const jacobians = getModel()->getJacobians();
  computeEigenValuesVectors(rightEv,leftEv,eValues,normal);

  m1 = jacobPlus + jacobMin;

  // CFout << "K = " << "\n" << m1 << "\n";
  CFout << "eValues = " << eValues << "\n";

  m2 = m1*rightEv;
  m3 = leftEv*m2;

  CFout << "LKR = " << "\n" << m3 << "\n";

  for (CFuint i = 0; i < 4; ++i) {
    for (CFuint j = 0; j < 4; ++j) {
      if (i != j) {
	if (std::abs(m3(i,j)) > 0.00000001) {
	  CFout << "eValues = " << eValues << "\n";
	  CFout << "LKR = " << "\n" << m3 << "\n";
	  throw Common::BadValueException (FromHere(),);
	}
      }
    }
  }

  m1 = (*jacobians)[0]*nCsi;
  m1 += (*jacobians)[1]*nEta;

  // CFout << "J = " << "\n" << m1 << "\n";

  m2 = m1*rightEv;
  m3 = leftEv*m2;

  CFout << "LJR = " << "\n" << m3 << "\n";

#endif
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DChar::computeScalarJacobian(const RealVector& normal,
				 RealVector& jacob)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal cosD = linearData[EulerTerm::VX]/linearData[EulerTerm::V];
  const CFreal sinD = linearData[EulerTerm::VY]/linearData[EulerTerm::V];
  const CFreal nCsi = normal[XX]*cosD + normal[YY]*sinD;
  const CFreal nEta = -normal[XX]*sinD + normal[YY]*cosD;
  const CFreal speed = linearData[EulerTerm::V];
  const CFreal M     = speed/linearData[EulerTerm::A];
  const CFreal M2    = M*M;
  const CFreal eps   = 0.05;
  const CFreal beta   = sqrt(max<CFreal>(eps*eps,std::abs(M2 - 1.)));
  const CFreal beta2  = beta*beta;
  const CFreal chi    = beta/max<CFreal>(M,1.);
  const CFreal nuPlus = 0.5*(M2 - 1. + beta2)/beta2;
  const CFreal um = speed*nCsi;

  if (jacob.size() == 2) {
    jacob[0] = um;
    jacob[1] = um;
  }
  else {
    // you get here if you are solving a fully supersonic flow
    // (four scalar equations)
    cf_assert(jacob.size() == 4);
    cf_assert(M >= 1.0);
    jacob[0] = speed*chi*(nuPlus*nCsi + nEta/beta);
    jacob[1] = speed*chi*(nuPlus*nCsi - nEta/beta);
    jacob[2] = um;
    jacob[3] = um;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
