#include "NavierStokes/NavierStokes.hh"
#include <numeric>

#include "Euler2DPrim.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DPrim, ConvectiveVarSet, NavierStokesModule, 1> euler2DPrimProvider("Euler2DPrim");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrim::Euler2DPrim(Common::SafePtr<BaseTerm> term) :
  Euler2DVarSet(term)
{
   vector<std::string> names(4);
   names[0] = "rho";
   names[1] = "u";
   names[2] = "v";
   names[3] = (!getModel()->isIncompressible()) ? "p" : "dp";
   setVarNames(names);
}
      
//////////////////////////////////////////////////////////////////////////////

Euler2DPrim::~Euler2DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrim::computeJacobians()
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avA   = linearData[EulerTerm::A];

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](0,0) = avU;
  (*jacobians)[0](0,1) = avRho;
  (*jacobians)[0](1,1) = avU;
  (*jacobians)[0](1,3) = 1.0/avRho;
  (*jacobians)[0](2,2) = avU;
  (*jacobians)[0](3,1) = avRho*avA*avA;
  (*jacobians)[0](3,3) = avU;

  (*jacobians)[1](0,0) = avV;
  (*jacobians)[1](0,2) = avRho;
  (*jacobians)[1](1,1) = avV;
  (*jacobians)[1](2,2) = avV;
  (*jacobians)[1](2,3) = 1.0/avRho;
  (*jacobians)[1](3,2) = avRho*avA*avA;
  (*jacobians)[1](3,3) = avV;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
          RealMatrix& leftEv,
          RealVector& eValues,
          const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal um = avU*nx+avV*ny;
  const CFreal rhoA = avRho*avA;

  rightEv(0,0) = 1.;
  rightEv(0,1) = 0.0;
  rightEv(0,2) = avRho/avA;
  rightEv(0,3) = avRho/avA;
  rightEv(1,0) = 0.0;
  rightEv(1,1) = -ny;
  rightEv(1,2) = nx;
  rightEv(1,3) = -nx;
  rightEv(2,0) = 0.0;
  rightEv(2,1) = nx;;
  rightEv(2,2) = ny;
  rightEv(2,3) = -ny;
  rightEv(3,0) = 0.0;
  rightEv(3,1) = 0.0;
  rightEv(3,2) = rhoA;
  rightEv(3,3) = rhoA;

  leftEv(0,0) = 1.;
  leftEv(0,1) = 0.0;
  leftEv(0,2) = 0.0;
  leftEv(0,3) = -1.0/(avA*avA);
  leftEv(1,0) = 0.0;
  leftEv(1,1) = -ny;
  leftEv(1,2) = nx;
  leftEv(1,3) = 0.0;
  leftEv(2,0) = 0.0;
  leftEv(2,1) = 0.5*nx;
  leftEv(2,2) = 0.5*ny;
  leftEv(2,3) = 0.5/rhoA;
  leftEv(3,0) = 0.0;
  leftEv(3,1) = -0.5*nx;
  leftEv(3,2) = -0.5*ny;
  leftEv(3,3) = 0.5/rhoA;

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

CFuint Euler2DPrim::getBlockSeparator() const
{
  return 4;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrim::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal Vn  = avU*nx + avV*ny;
  const CFreal a   = linearData[EulerTerm::A];

  eValues[0] = Vn;
  eValues[1] = Vn;
  eValues[2] = Vn + a;
  eValues[3] = Vn - a;
  for (CFuint iEq = 0; iEq < 4; ++iEq) {
    _eValuesP[iEq] = max(0.,eValues[iEq]);
    _eValuesM[iEq] = min(0.,eValues[iEq]);
  }

  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal a2 = a*a;
  const CFreal rhoa = rho*a;
  const CFreal rhoOverA = rho/a;
  const CFreal nxny = nx*ny;
  const CFreal nx2 = nx*nx;
  const CFreal ny2 = ny*ny;

  CFreal e1 = _eValuesP[0];
  CFreal e2 = _eValuesP[1];
  CFreal e3Min4 = 0.5*(_eValuesP[2] - _eValuesP[3]);
  CFreal e3Plus4 = 0.5*(_eValuesP[2] + _eValuesP[3]);
  CFreal e3Min4nx = e3Min4*nx;
  CFreal e3Min4ny = e3Min4*ny;

  jacobPlus(0,0) = e1;
  jacobPlus(0,1) = e3Min4nx*rhoOverA;
  jacobPlus(0,2) = e3Min4ny*rhoOverA;
  jacobPlus(0,3) = (e3Plus4 - e1)/a2;
  jacobPlus(1,1) = e2*ny2 + e3Plus4*nx2;
  jacobPlus(1,2) = -e2*nxny + e3Plus4*nxny;
  jacobPlus(1,3) = e3Min4nx/rhoa;
  jacobPlus(2,1) = -e2*nxny + e3Plus4*nxny;
  jacobPlus(2,2) = e2*nx2 + e3Plus4*ny2;
  jacobPlus(2,3) = e3Min4ny/rhoa;
  jacobPlus(3,1) = e3Min4nx*rhoa;
  jacobPlus(3,2) = e3Min4ny*rhoa;
  jacobPlus(3,3) = e3Plus4;

  e1 = _eValuesM[0];
  e2 = _eValuesM[1];
  e3Min4 = 0.5*(_eValuesM[2] - _eValuesM[3]);
  e3Plus4 = 0.5*(_eValuesM[2] + _eValuesM[3]);
  e3Min4nx = e3Min4*nx;
  e3Min4ny = e3Min4*ny;

  jacobMin(0,0) = e1;
  jacobMin(0,1) = e3Min4nx*rhoOverA;
  jacobMin(0,2) = e3Min4ny*rhoOverA;
  jacobMin(0,3) = (e3Plus4 - e1)/a2;
  jacobMin(1,1) = e2*ny2 + e3Plus4*nx2;
  jacobMin(1,2) = -e2*nxny + e3Plus4*nxny;
  jacobMin(1,3) = e3Min4nx/rhoa;
  jacobMin(2,1) = -e2*nxny + e3Plus4*nxny;
  jacobMin(2,2) = e2*nx2 + e3Plus4*ny2;
  jacobMin(2,3) = e3Min4ny/rhoa;
  jacobMin(3,1) = e3Min4nx*rhoa;
  jacobMin(3,2) = e3Min4ny*rhoa;
  jacobMin(3,3) = e3Plus4;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrim::computePhysicalData(const State& state, RealVector& data)
  
{
  const CFreal rho = state[0]; // this model cannot treat the pure incompresible case
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal p = getModel()->getPressureFromState(state[3]);
  const CFreal V2 = u*u + v*v;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal pOvRho = p/rho;

  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = state[3];
  data[EulerTerm::H] = gammaDivGammaMinus1*pOvRho + 0.5*V2;
  data[EulerTerm::E] = data[EulerTerm::H] - pOvRho;
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = pOvRho/getModel()->getR();
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::GAMMA] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrim::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
  state[0] = data[EulerTerm::RHO];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DPrim::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrim::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::RHO];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrim::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::RHO];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/refData[EulerTerm::P];
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
