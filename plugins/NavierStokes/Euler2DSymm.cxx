#include "NavierStokes/NavierStokes.hh"
#include "Euler2DSymm.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DSymm, ConvectiveVarSet, NavierStokesModule, 1>
euler2DSymmProvider("Euler2DSymm");

//////////////////////////////////////////////////////////////////////////////

Euler2DSymm::Euler2DSymm(Common::SafePtr<BaseTerm> term) :
  Euler2DVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DSymm::~Euler2DSymm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSymm::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal avA = linearData[EulerTerm::A];
  const CFreal speed = sqrt(avU*avU + avV+avV);

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](0,0) = speed;
  (*jacobians)[0](0,1) = avA;
  (*jacobians)[0](1,0) = avA;
  (*jacobians)[0](1,1) = speed;
  (*jacobians)[0](2,2) = speed;
  (*jacobians)[0](3,3) = speed;

  (*jacobians)[1](0,2) = avA;
  (*jacobians)[1](2,0) = avA;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSymm::computeEigenValuesVectors(RealMatrix& rightEv,
					RealMatrix& leftEv,
					RealVector& eValues,
					const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal speed = sqrt(avU*avU + avV*avV); //absolute velocity in the dof !!!
  cf_assert(speed > 0.);
  const CFreal cosD = avU/speed;
  const CFreal sinD = avV/speed;
  const CFreal nCsi = nx*cosD + ny*sinD;
  const CFreal nEta = -nx*sinD + ny*cosD;
  const CFreal um = avU*nx + avV*ny;
  //  const CFreal um = speed*nCsi;

  //the eigenvalues of the matrix Aq*nCsi + Bq*nEta are the same as the ones
  //for  A*nx + B*ny
  if (rightEv.nbRows() == 3) {
    rightEv(0,1) = 1.0;
    rightEv(0,2) = 1.0;
    rightEv(1,0) = -nEta;
    rightEv(1,1) = nCsi;
    rightEv(1,2) = -nCsi;
    rightEv(2,0) = nCsi;
    rightEv(2,1) = nEta;
    rightEv(2,2) = -nEta;

    leftEv(0,1) = -nEta;
    leftEv(0,2) = nCsi;
    leftEv(1,0) = 0.5;
    leftEv(1,1) = 0.5*nCsi;
    leftEv(1,2) = 0.5*nEta;
    leftEv(2,0) = 0.5;
    leftEv(2,1) = -0.5*nCsi;
    leftEv(2,2) = -0.5*nEta;

    eValues[0] = um;
    eValues[1] = um + avA;
    eValues[2] = um - avA;
  }
  else {
    cf_assert(rightEv.nbRows() == 4);
    
    rightEv(0,2) = 1.0;
    rightEv(0,3) = 1.0;
    rightEv(1,0) = -nEta;
    rightEv(1,1) = -nEta;
    rightEv(1,2) = nCsi;
    rightEv(1,3) = -nCsi;
    rightEv(2,0) = nCsi;
    rightEv(2,1) = nCsi;
    rightEv(2,2) = nEta;
    rightEv(2,3) = -nEta;
    rightEv(3,0) = 1.0;

    leftEv(0,3) = 1.0;
    leftEv(1,1) = -nEta;
    leftEv(1,2) = nCsi;
    leftEv(1,3) = -1.0;
    leftEv(2,0) = 0.5;
    leftEv(2,1) = 0.5*nCsi;
    leftEv(2,2) = 0.5*nEta;
    leftEv(3,0) = 0.5;
    leftEv(3,1) = -0.5*nCsi;
    leftEv(3,2) = -0.5*nEta;

    eValues[0] = um;
    eValues[1] = um;
    eValues[2] = um + avA;
    eValues[3] = um - avA;
  }

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
	<< rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
	<< leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DChar::computeEigenValuesVectors" << "\n"
	<< eValues << "\n" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DSymm::getBlockSeparator() const
{
  return 3;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSymm::splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal speed = sqrt(avU*avU + avV*avV);
  const CFreal cosD = avU/speed;
  const CFreal sinD = avV/speed;
  const CFreal nCsi = nx*cosD + ny*sinD;
  const CFreal nEta = -nx*sinD + ny*cosD;
  const CFreal nEtaCsi = nEta*nCsi;
  const CFreal nCsi2 = nCsi*nCsi;
  const CFreal nEta2 = nEta*nEta;
  const CFreal Vn  = avU*nx + avV*ny;
  const CFuint sizeJacob = jacobPlus.nbRows();

  cf_assert(eValues.size() == 4);
  eValues[0] = eValues[1] = Vn;
  eValues[2] = Vn + avA;
  eValues[3] = Vn - avA;

  for (CFuint iEq = 0; iEq < 4; ++iEq) {
    _eValuesP[iEq] = max(0.,eValues[iEq]);
    _eValuesM[iEq] = min(0.,eValues[iEq]);
  }

  CFreal e2 = _eValuesP[1];
  CFreal e3Min4 = 0.5*(_eValuesP[2] - _eValuesP[3]);
  CFreal e3Plus4 = 0.5*(_eValuesP[2] + _eValuesP[3]);
  CFreal e3Min4nEta = e3Min4*nEta;
  CFreal e3Min4nCsi = e3Min4*nCsi;
  CFreal e234 = -e2 + e3Plus4;

  if (sizeJacob == 3) {
    jacobPlus(0,0) = e3Plus4;
    jacobPlus(0,1) = e3Min4*nCsi;
    jacobPlus(0,2) = e3Min4*nEta;
    jacobPlus(1,0) = e3Min4*nCsi;
    jacobPlus(1,1) = e2*nEta2 + e3Plus4*nCsi2;
    jacobPlus(1,2) = nEtaCsi*e234;
    jacobPlus(2,0) = e3Min4*nEta;
    jacobPlus(2,1) = nEtaCsi*e234;
    jacobPlus(2,2) = e2*nCsi2 + e3Plus4*nEta2;

    //  jacobPlus[0] = e3Plus4;
    //  jacobPlus[1] = e3Min4*nCsi;
    //  jacobPlus[2] = e3Min4*nEta;

    //  jacobPlus[4] = e3Min4*nCsi;
    //  jacobPlus[5] = e2*nEta2 + e3Plus4*nCsi2;
    //  jacobPlus[6] = nEtaCsi*e234;

    //  jacobPlus[8] = e3Min4*nEta;
    //  jacobPlus[9] = nEtaCsi*e234;
    //  jacobPlus[10] = e2*nCsi2 + e3Plus4*nEta2;

    e2 = _eValuesM[1];
    e3Min4 = 0.5*(_eValuesM[2] - _eValuesM[3]);
    e3Plus4 = 0.5*(_eValuesM[2] + _eValuesM[3]);
    e3Min4nEta = e3Min4*nEta;
    e3Min4nCsi = e3Min4*nCsi;
    e234 = -e2 + e3Plus4;

    jacobMin(0,0) = e3Plus4;
    jacobMin(0,1) = e3Min4*nCsi;
    jacobMin(0,2) = e3Min4*nEta;
    jacobMin(1,0) = e3Min4*nCsi;
    jacobMin(1,1) = e2*nEta2 + e3Plus4*nCsi2;
    jacobMin(1,2) = nEtaCsi*e234;
    jacobMin(2,0) = e3Min4*nEta;
    jacobMin(2,1) = nEtaCsi*e234;
    jacobMin(2,2) = e2*nCsi2 + e3Plus4*nEta2;
  }
  else {
    cf_assert(sizeJacob == 4);

    jacobPlus(0,0) = e3Plus4;
    jacobPlus(0,1) = e3Min4*nCsi;
    jacobPlus(0,2) = e3Min4*nEta;
    jacobPlus(1,0) = e3Min4*nCsi;
    jacobPlus(1,1) = e2*nEta2 + e3Plus4*nCsi2;
    jacobPlus(1,2) = nEtaCsi*e234;
    jacobPlus(2,0) = e3Min4*nEta;
    jacobPlus(2,1) = nEtaCsi*e234;
    jacobPlus(2,2) = e2*nCsi2 + e3Plus4*nEta2;
    jacobPlus(3,3) = e2; //e2 = e1

    e2 = _eValuesM[1];
    e3Min4 = 0.5*(_eValuesM[2] - _eValuesM[3]);
    e3Plus4 = 0.5*(_eValuesM[2] + _eValuesM[3]);
    e3Min4nEta = e3Min4*nEta;
    e3Min4nCsi = e3Min4*nCsi;
    e234 = -e2 + e3Plus4;

    jacobMin(0,0) = e3Plus4;
    jacobMin(0,1) = e3Min4*nCsi;
    jacobMin(0,2) = e3Min4*nEta;
    jacobMin(1,0) = e3Min4*nCsi;
    jacobMin(1,1) = e2*nEta2 + e3Plus4*nCsi2;
    jacobMin(1,2) = nEtaCsi*e234;
    jacobMin(2,0) = e3Min4*nEta;
    jacobMin(2,1) = nEtaCsi*e234;
    jacobMin(2,2) = e2*nCsi2 + e3Plus4*nEta2;
    jacobMin(3,3) = e2; //e2 = e1
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSymm::computeScalarJacobian(const RealVector& normal,
					RealVector& jacob)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  
  jacob[0] = linearData[EulerTerm::VX]*normal[XX] +
    linearData[EulerTerm::VY]*normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
