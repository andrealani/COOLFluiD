#include "NavierStokes/NavierStokes.hh"
#include "Euler3DCons.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DCons, ConvectiveVarSet, NavierStokesModule, 1>
euler3DConsProvider("Euler3DCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DCons::Euler3DCons(Common::SafePtr<BaseTerm> term) :
  Euler3DVarSet(term),
  _rightEv(5,5),
  _leftEv(5,5)
{

  vector<std::string> names(5);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoW";
  names[4] = "rhoE";
  setVarNames(names);
}


//////////////////////////////////////////////////////////////////////////////

Euler3DCons::~Euler3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setup()
{
  Euler3DVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setConstJacob()
{
  vector<RealMatrix>* const jacobians = PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  (*jacobians)[0](0,1) = 1.0;
  (*jacobians)[0](1,4) = gammaMinus1;

  (*jacobians)[1](0,2) = 1.0;
  (*jacobians)[1](2,4) = gammaMinus1;

  (*jacobians)[2](0,3) = 1.0;
  (*jacobians)[2](3,4) = gammaMinus1;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal avW = linearData[EulerTerm::VZ];
  const CFreal avH = linearData[EulerTerm::H];
  const CFreal uu = avU*avU;
  const CFreal vv = avV*avV;
  const CFreal ww = avW*avW;
  const CFreal uv = avU*avV;
  const CFreal uw = avU*avW;
  const CFreal vw = avV*avW;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus3 = gamma - 3.;
  const CFreal sumVel2 = 0.5*gammaMinus1*(uu + vv + ww);

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](1,0) = sumVel2 - uu;
  (*jacobians)[0](1,1) = -gammaMinus3*avU;
  (*jacobians)[0](1,2) = -gammaMinus1*avV;
  (*jacobians)[0](1,3) = -gammaMinus1*avW;
  (*jacobians)[0](2,0) = -uv;
  (*jacobians)[0](2,1) = avV;
  (*jacobians)[0](2,2) = avU;
  (*jacobians)[0](3,0) = -uw;
  (*jacobians)[0](3,1) = avW;
  (*jacobians)[0](3,3) = avU;
  (*jacobians)[0](4,0) = avU*(sumVel2 - avH);
  (*jacobians)[0](4,1) = -gammaMinus1*uu + avH;
  (*jacobians)[0](4,2) = -gammaMinus1*uv;
  (*jacobians)[0](4,3) = -gammaMinus1*uw;
  (*jacobians)[0](4,4) = gamma*avU;

  (*jacobians)[1](1,0) = -uv;
  (*jacobians)[1](1,1) = avV;
  (*jacobians)[1](1,2) = avU;
  (*jacobians)[1](2,0) = sumVel2 - vv;
  (*jacobians)[1](2,1) = -gammaMinus1*avU;
  (*jacobians)[1](2,2) = -gammaMinus3*avV;
  (*jacobians)[1](2,3) = -gammaMinus1*avW;
  (*jacobians)[1](3,0) = -vw;
  (*jacobians)[1](3,2) = avW;
  (*jacobians)[1](3,3) = avV;
  (*jacobians)[1](4,0) = avV*(sumVel2 - avH);
  (*jacobians)[1](4,1) = -gammaMinus1*uv;
  (*jacobians)[1](4,2) = -gammaMinus1*vv + avH;
  (*jacobians)[1](4,3) = -gammaMinus1*vw;
  (*jacobians)[1](4,4) = gamma*avV;

  (*jacobians)[2](1,0) = -uw;
  (*jacobians)[2](1,1) = avW;
  (*jacobians)[2](1,3) = avU;
  (*jacobians)[2](2,0) = -vw;
  (*jacobians)[2](2,2) = avW;
  (*jacobians)[2](2,3) = avV;
  (*jacobians)[2](3,0) = sumVel2 - ww;
  (*jacobians)[2](3,1) = -gammaMinus1*avU;
  (*jacobians)[2](3,2) = -gammaMinus1*avV;
  (*jacobians)[2](3,3) = -gammaMinus3*avW;
  (*jacobians)[2](4,0) = avW*(sumVel2 - avH);
  (*jacobians)[2](4,1) = -gammaMinus1*uw;
  (*jacobians)[2](4,2) = -gammaMinus1*vw;
  (*jacobians)[2](4,3) = -gammaMinus1*ww + avH;
  (*jacobians)[2](4,4) = gamma*avW;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                        RealMatrix& leftEv,
                                        RealVector& eValues,
                                        const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avW   = linearData[EulerTerm::VZ];
  const CFreal avH   = linearData[EulerTerm::H];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal avK   = 0.5*linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFreal um = avU*nx + avV*ny + avW*nz;
  const CFreal ra = 0.5*avRho/avA;
  const CFreal avA2 = avA*avA;
  const CFreal invAvRho = 1./avRho;
  const CFreal k1 = gammaMinus1*avK/avA2;
  const CFreal k2 = 1.0 - k1;
  const CFreal k3 = -gammaMinus1/avA2;
  const CFreal uDivA = gammaMinus1*avU/avA;
  const CFreal vDivA = gammaMinus1*avV/avA;
  const CFreal wDivA = gammaMinus1*avW/avA;
  const CFreal uDivA2 = uDivA/avA;
  const CFreal vDivA2 = vDivA/avA;
  const CFreal wDivA2 = wDivA/avA;
  const CFreal rhoA = avRho*avA;

  rightEv(0,0) = nx;
  rightEv(0,1) = ny;
  rightEv(0,2) = nz;
  rightEv(0,3) = ra;
  rightEv(0,4) = ra;
  rightEv(1,0) = avU*nx;
  rightEv(1,1) = avU*ny - avRho*nz;
  rightEv(1,2) = avU*nz + avRho*ny;
  rightEv(1,3) = ra*(avU + avA*nx);
  rightEv(1,4) = ra*(avU - avA*nx);
  rightEv(2,0) = avV*nx + avRho*nz;
  rightEv(2,1) = avV*ny;
  rightEv(2,2) = avV*nz - avRho*nx;
  rightEv(2,3) = ra*(avV + avA*ny);
  rightEv(2,4) = ra*(avV - avA*ny);
  rightEv(3,0) = avW*nx - avRho*ny;
  rightEv(3,1) = avW*ny + avRho*nx;
  rightEv(3,2) = avW*nz;
  rightEv(3,3) = ra*(avW + avA*nz);
  rightEv(3,4) = ra*(avW - avA*nz);
  rightEv(4,0) = avK*nx + avRho*(avV*nz - avW*ny);
  rightEv(4,1) = avK*ny + avRho*(avW*nx - avU*nz);
  rightEv(4,2) = avK*nz + avRho*(avU*ny - avV*nx);
  rightEv(4,3) = ra*(avH + avA*um);
  rightEv(4,4) = ra*(avH - avA*um);

  leftEv(0,0) = nx*k2 - invAvRho*(avV*nz - avW*ny);
  leftEv(0,1) = uDivA2*nx;
  leftEv(0,2) = vDivA2*nx + nz*invAvRho;
  leftEv(0,3) = wDivA2*nx - ny*invAvRho;
  leftEv(0,4) = k3*nx;
  leftEv(1,0) = ny*k2 - invAvRho*(avW*nx - avU*nz);
  leftEv(1,1) = uDivA2*ny - nz*invAvRho;
  leftEv(1,2) = vDivA2*ny;
  leftEv(1,3) = wDivA2*ny + nx*invAvRho;
  leftEv(1,4) = k3*ny;
  leftEv(2,0) = nz*k2 - invAvRho*(avU*ny - avV*nx);
  leftEv(2,1) = uDivA2*nz + ny*invAvRho;
  leftEv(2,2) = vDivA2*nz - nx*invAvRho;
  leftEv(2,3) = wDivA2*nz;
  leftEv(2,4) = k3*nz;
  leftEv(3,0) = avA*invAvRho*(k1 - um/avA);
  leftEv(3,1) = invAvRho*(nx - uDivA);
  leftEv(3,2) = invAvRho*(ny - vDivA);
  leftEv(3,3) = invAvRho*(nz - wDivA);
  leftEv(3,4) = gammaMinus1/rhoA;
  leftEv(4,0) = avA*invAvRho*(k1 + um/avA);
  leftEv(4,1) = invAvRho*(-nx - uDivA);
  leftEv(4,2) = invAvRho*(-ny - vDivA);
  leftEv(4,3) = invAvRho*(-nz - wDivA);
  leftEv(4,4) = gammaMinus1/rhoA;

  eValues[0] = um;
  eValues[1] = um;
  eValues[2] = um;
  eValues[3] = um + avA;
  eValues[4] = um - avA;

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler3DCons::computeEigenValuesVectors" << "\n"
        << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler3DCons::computeEigenValuesVectors" << "\n"
        << leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler3DCons::computeEigenValuesVectors" << "\n"
		 << eValues << "\n" << "\n");


 //  // test !!!!!
//   computeJacobians();
//   vector<RealMatrix>* const jacobians =
//     PhysicalModelStack::getActive()->getImplementor()->getJacobians();

//   RealMatrix mat1(5,5);
//   RealMatrix mat2(5,5);
//   mat1 = nx*(*jacobians)[0] + ny*(*jacobians)[1] + nz*(*jacobians)[2];
//   mat2 = leftEv*(mat1*rightEv);

//   //  mat = leftEv*(((*jacobians)[0]*nx + (*jacobians)[1]*ny + (*jacobians)[2]*nz)*rightEv);
//   CFLogInfo( "The Matrix: \n" << mat2 << "\n\n");

//   CFLogInfo("EigenValues = " << eValues << "\n\n\n");

//   //   CFLog(3, "\n" << "Identity: " << "\n"
//   //           << leftEv*rightEv << "\n"
//   //           << "Identity: " << "\n" <<
//   //           rightEv*leftEv << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::splitJacobian(RealMatrix& jacobPlus,
                             RealMatrix& jacobMin,
                             RealVector& eValues,
                             const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal nx = normal[0];
  const CFreal ny = normal[1];
  const CFreal nz = normal[2];
  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avW   = linearData[EulerTerm::VZ];
  const CFreal avH   = linearData[EulerTerm::H];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal avK   = 0.5*linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFreal um = avU*nx + avV*ny + avW*nz;
  const CFreal ra = 0.5*avRho/avA;
  const CFreal avA2 = avA*avA;
  const CFreal invAvRho = 1./avRho;
  const CFreal k1 = gammaMinus1*avK/avA2;
  const CFreal k2 = 1.0 - k1;
  const CFreal k3 = -gammaMinus1/avA2;
  const CFreal uDivA = gammaMinus1*avU/avA;
  const CFreal vDivA = gammaMinus1*avV/avA;
  const CFreal wDivA = gammaMinus1*avW/avA;
  const CFreal uDivA2 = uDivA/avA;
  const CFreal vDivA2 = vDivA/avA;
  const CFreal wDivA2 = wDivA/avA;
  const CFreal rhoA = avRho*avA;

  _rightEv(0,0) = nx;
  _rightEv(0,1) = ny;
  _rightEv(0,2) = nz;
  _rightEv(0,3) = ra;
  _rightEv(0,4) = ra;

  _rightEv(1,0) = avU*nx;
  _rightEv(1,1) = avU*ny - avRho*nz;
  _rightEv(1,2) = avU*nz + avRho*ny;
  _rightEv(1,3) = ra*(avU + avA*nx);
  _rightEv(1,4) = ra*(avU - avA*nx);

  _rightEv(2,0) = avV*nx + avRho*nz;
  _rightEv(2,1) = avV*ny;
  _rightEv(2,2) = avV*nz - avRho*nx;
  _rightEv(2,3) = ra*(avV + avA*ny);
  _rightEv(2,4) = ra*(avV - avA*ny);

  _rightEv(3,0) = avW*nx - avRho*ny;
  _rightEv(3,1) = avW*ny + avRho*nx;
  _rightEv(3,2) = avW*nz;
  _rightEv(3,3) = ra*(avW + avA*nz);
  _rightEv(3,4) = ra*(avW - avA*nz);

  _rightEv(4,0) = avK*nx + avRho*(avV*nz - avW*ny);
  _rightEv(4,1) = avK*ny + avRho*(avW*nx - avU*nz);
  _rightEv(4,2) = avK*nz + avRho*(avU*ny - avV*nx);
  _rightEv(4,3) = ra*(avH + avA*um);
  _rightEv(4,4) = ra*(avH - avA*um);

  _leftEv(0,0) = nx*k2 - invAvRho*(avV*nz - avW*ny);
  _leftEv(0,1) = uDivA2*nx;
  _leftEv(0,2) = vDivA2*nx + nz*invAvRho;
  _leftEv(0,3) = wDivA2*nx - ny*invAvRho;
  _leftEv(0,4) = k3*nx;
  _leftEv(1,0) = ny*k2 - invAvRho*(avW*nx - avU*nz);
  _leftEv(1,1) = uDivA2*ny - nz*invAvRho;
  _leftEv(1,2) = vDivA2*ny;
  _leftEv(1,3) = wDivA2*ny + nx*invAvRho;
  _leftEv(1,4) = k3*ny;
  _leftEv(2,0) = nz*k2 - invAvRho*(avU*ny - avV*nx);
  _leftEv(2,1) = uDivA2*nz + ny*invAvRho;
  _leftEv(2,2) = vDivA2*nz - nx*invAvRho;
  _leftEv(2,3) = wDivA2*nz;
  _leftEv(2,4) = k3*nz;
  _leftEv(3,0) = avA*invAvRho*(k1 - um/avA);
  _leftEv(3,1) = invAvRho*(nx - uDivA);
  _leftEv(3,2) = invAvRho*(ny - vDivA);
  _leftEv(3,3) = invAvRho*(nz - wDivA);
  _leftEv(3,4) = gammaMinus1/rhoA;
  _leftEv(4,0) = avA*invAvRho*(k1 + um/avA);
  _leftEv(4,1) = invAvRho*(-nx - uDivA);
  _leftEv(4,2) = invAvRho*(-ny - vDivA);
  _leftEv(4,3) = invAvRho*(-nz - wDivA);
  _leftEv(4,4) = gammaMinus1/rhoA;

  eValues[0] = um;
  eValues[1] = um;
  eValues[2] = um;
  eValues[3] = um + avA;
  eValues[4] = um - avA;

  // compute the eigen values + and -
  for (CFuint iEq = 0; iEq < 5; ++iEq) {
    _eValuesP[iEq] = max(0.,eValues[iEq]);
    _eValuesM[iEq] = min(0.,eValues[iEq]);
  }

  // compute jacobian + and -
  jacobPlus = _rightEv*(_eValuesP*_leftEv);
  jacobMin  = _rightEv*(_eValuesM*_leftEv);

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler3DCons::computeEigenValuesVectors" << "\n"
        << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler3DCons::computeEigenValuesVectors" << "\n"
        << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler3DCons::computeEigenValuesVectors" << "\n"
        << eValues << "\n" << "\n");


  // test !!!!!
  // computeJacobians();
  //   vector<RealMatrix>* const jacobians =
  //     PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  //   RealMatrix mat1(5,5);
  //   RealMatrix mat2(5,5);
  //   RealMatrix mat3(5,5);
  //   mat1 = nx*(*jacobians)[0] + ny*(*jacobians)[1] + nz*(*jacobians)[2];
  //   mat2 = mat1*_rightEv;
  //   mat3 = _leftEv*mat2;

  //   CFLogInfo( "The Matrix: \n" << mat3 << "\n\n");

  //   CFLogInfo("EigenValues = " << eValues << "\n\n\n");

  //   mat1 = _leftEv*_rightEv;
  //   cout << "ID1 = \n"<< mat1 << endl;
  //   mat1 = _rightEv*_leftEv;
  //   cout << "ID2 = \n"<< mat1 << endl;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setEigenVect1(RealVector& r1,
                                State& state,
                                const RealVector& normal)
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal w = state[3]/state[0];
  const CFreal nx = normal[0];
  const CFreal ny = normal[1];
  const CFreal nz = normal[2];

  r1[0] = nx;
  r1[1] = u*nx;
  r1[2] = v*nx + state[0]*nz;
  r1[3] = w*nx - state[0]*ny;
  r1[4] = 0.5*(u*u + v*v + w*w) + state[0]*(v*nz - w*ny);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setEigenVect2(RealVector& r2,
                                State& state,
                                const RealVector& normal)
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal w = state[3]/state[0];
  const CFreal nx = normal[0];
  const CFreal ny = normal[1];
  const CFreal nz = normal[2];

  r2[0] = ny;
  r2[1] = u*ny - state[0]*nz;
  r2[2] = v*ny;
  r2[3] = w*ny + state[0]*nx;
  r2[4] = 0.5*(u*u + v*v + w*w) + state[0]*(w*nx - u*nz);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setEigenVect3(RealVector& r3,
                                State& state,
                                const RealVector& normal)
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal w = state[3]/state[0];
  const CFreal nx = normal[0];
  const CFreal ny = normal[1];
  const CFreal nz = normal[2];

  r3[0] = nz;
  r3[1] = u*nz + state[0]*ny;
  r3[2] = v*nz - state[0]*nx;
  r3[3] = w*nz;
  r3[4] = 0.5*(u*u + v*v + w*w) + state[0]*(u*ny - v*nx);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setEigenVect4(RealVector& r4,
                                State& state,
                                const RealVector& normal)
{
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = 1./(gamma - 1.);
  const CFreal rhoInv = 1./state[0];
  const CFreal u = state[1]*rhoInv;
  const CFreal v = state[2]*rhoInv;
  const CFreal w = state[3]*rhoInv;
  const CFreal nx = normal[0];
  const CFreal ny = normal[1];
  const CFreal nz = normal[2];
  const CFreal k  =  0.5*(u*u + v*v + w*w);
  const CFreal p = gammaMinus1*(state[4] - state[0]*k);
  const CFreal a = sqrt(gamma*p*rhoInv);
  const CFreal H = rhoInv*(gammaDivGammaMinus1*p + state[0]*k);
  const CFreal ra = 0.5*state[0]/a;

  r4[0] = ra;
  r4[1] = ra*(u + a*nx);
  r4[2] = ra*(v + a*ny);
  r4[3] = ra*(w + a*nz);
  r4[4] = ra*(H + a*(u*nx + v*ny + w*nz));
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setEigenVect5(RealVector& r5,
                                State& state,
                                const RealVector& normal)
{
  const CFreal gamma = getModel()->getGamma();
 const CFreal gammaMinus1 = gamma - 1.;
 const CFreal gammaDivGammaMinus1 = 1./(gamma - 1.);
  const CFreal rhoInv = 1./state[0];
  const CFreal u = state[1]*rhoInv;
  const CFreal v = state[2]*rhoInv;
  const CFreal w = state[3]*rhoInv;
  const CFreal nx = normal[0];
  const CFreal ny = normal[1];
  const CFreal nz = normal[2];
  const CFreal k  =  0.5*(u*u + v*v + w*w);
  const CFreal p = gammaMinus1*(state[4] - state[0]*k);
  const CFreal a = sqrt(gamma*p*rhoInv);
  const CFreal H = rhoInv*(gammaDivGammaMinus1*p + state[0]*k);
  const CFreal ra = 0.5*state[0]/a;

  r5[0] = ra;
  r5[1] = ra*(u - a*nx);
  r5[2] = ra*(v - a*ny);
  r5[3] = ra*(w - a*nz);
  r5[4] = ra*(H - a*(u*nx + v*ny + w*nz));
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::computePhysicalData(const State& state, RealVector& data)
{
  // we assume that if conservative variables are used, the flow is compressible
  // p = static pressure
  const CFreal rho  = state[0];
  const CFreal ovRho = 1./rho;
  const CFreal u = state[1]*ovRho;
  const CFreal v = state[2]*ovRho;
  const CFreal w = state[3]*ovRho;
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = gammaMinus1*(state[4] - 0.5*rho*V2);

  const CFreal pOvRho = data[EulerTerm::P]*ovRho;
  data[EulerTerm::E] = state[4]*ovRho;
  data[EulerTerm::H] = data[EulerTerm::E] + pOvRho;
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = pOvRho/getModel()->getR();
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::VZ] = w;
  data[EulerTerm::GAMMA] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::computeStateFromPhysicalData(const RealVector& data,
					  State& state)
{
  const CFreal rho = data[EulerTerm::RHO];
  state[0] = rho;
  state[1] = rho*data[EulerTerm::VX];
  state[2] = rho*data[EulerTerm::VY];
  state[3] = rho*data[EulerTerm::VZ];
  state[4] = rho*data[EulerTerm::H] - data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DCons::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal w = state[3]/state[0];
  return sqrt(u*u + v*v + w*w);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData =
    getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::RHO];
  result[1] = state[1]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[4] = state[4]*refData[EulerTerm::RHO]*refData[EulerTerm::H];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData =
    getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::RHO];
  result[1] = state[1]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[2] = state[2]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[3] = state[3]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[4] = state[4]/(refData[EulerTerm::RHO]*refData[EulerTerm::H]);
}

//////////////////////////////////////////////////////////////////////////////
      
void Euler3DCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
  const RealVector& linearData = getModel()->getPhysicalData();
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal avW = linearData[EulerTerm::VZ];
  const CFreal avH = linearData[EulerTerm::H];
  const CFreal uu = avU*avU;
  const CFreal vv = avV*avV;
  const CFreal ww = avW*avW;
  const CFreal uv = avU*avV;
  const CFreal uw = avU*avW;
  const CFreal vw = avV*avW;
  const CFreal un = avU*nx + avV*ny + avW*nz;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus3 = gamma - 3.;
  const CFreal sumVel2 = 0.5*gammaMinus1*(uu + vv + ww);
  
  // constant entries 
  jacob(0,1) = nx;
  jacob(0,2) = ny;
  jacob(0,3) = nz;
  
  jacob(1,0) = (sumVel2 - uu)*nx -uv*ny -uw*nz;
  jacob(1,1) = -gammaMinus3*avU*nx + avV*ny + avW*nz;
  jacob(1,2) = -gammaMinus1*avV*nx + avU*ny;
  jacob(1,3) = -gammaMinus1*avW*nx + avU*nz;
  jacob(1,4) = gammaMinus1*nx;
  
  jacob(2,0) = -uv*nx + (sumVel2 - vv)*ny -vw*nz;
  jacob(2,1) = avV*nx -gammaMinus1*avU*ny;
  jacob(2,2) = avU*nx -gammaMinus3*avV*ny + avW*nz;
  jacob(2,3) = -gammaMinus1*avW*ny + avV*nz;
  jacob(2,4) = gammaMinus1*ny;
  
  jacob(3,0) = -uw*nx -vw*ny + (sumVel2 - ww)*nz;
  jacob(3,1) = avW*nx -gammaMinus1*avU*nz;
  jacob(3,2) = avW*ny -gammaMinus1*avV*nz;
  jacob(3,3) = avU*nx + avV*ny -gammaMinus3*avW*nz;
  jacob(3,4) = gammaMinus1*nz;
  
  jacob(4,0) = (sumVel2 - avH)*un;
  jacob(4,1) = (-gammaMinus1*uu + avH)*nx -gammaMinus1*uv*ny -gammaMinus1*uw*nz;
  jacob(4,2) = -gammaMinus1*uv*nx + (-gammaMinus1*vv + avH)*ny -gammaMinus1*vw*nz;
  jacob(4,3) = -gammaMinus1*uw*nx -gammaMinus1*vw*ny + (-gammaMinus1*ww + avH)*nz;
  jacob(4,4) = gamma*un;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
