#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPuvt.hh"
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

Environment::ObjectProvider<Euler2DPuvt, ConvectiveVarSet, NavierStokesModule, 1> 
euler2DPuvtProvider("Euler2DPuvt");
      
//////////////////////////////////////////////////////////////////////////////

Euler2DPuvt::Euler2DPuvt(Common::SafePtr<BaseTerm> term) :
  Euler2DVarSet(term)
{
  vector<std::string> names(4);
  names[0] = (!getModel()->isIncompressible()) ? "p" : "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvt::~Euler2DPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler2DPuvt::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DPuvt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DPuvt::getBlockSeparator() const
{
  return 4;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::splitJacobian(RealMatrix& jacobPlus,
				RealMatrix& jacobMin,
				RealVector& eValues,
				const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DPuvt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::computePhysicalData(const State& state, RealVector& data)
{
  const CFreal p = getModel()->getPressureFromState(state[0]);
  cf_assert(p > 0);
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal T = state[3];
  const CFreal rho = getModel()->getDensity(p,T);
  const CFreal V2 = u*u + v*v;
  const CFreal gamma = this->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.);
  const CFreal pOvRho = p/rho;
  
  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = state[0]; // p or dp
  data[EulerTerm::H] = gammaDivGammaMinus1*pOvRho + 0.5*V2;
  data[EulerTerm::E] = data[EulerTerm::H] - pOvRho;
  cf_assert(gamma*pOvRho > 0);
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = T;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::GAMMA] = gamma;
  
  if (this->_physDataNeedCoordinates) {
    const RealVector& node = state.getCoordinates();
    data[EulerTerm::XP] = node[XX];
    data[EulerTerm::YP] = node[YY];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::computeStateFromPhysicalData(const RealVector& data, State& state)
{
  state[0] = data[EulerTerm::P]; // p or dp
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DPuvt::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::P];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::setAdimensionalValues(const State& state,
                                        RealVector& result)
{  
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::P];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/(getModel()->getTempRef());
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvt::computePressureDerivatives(const Framework::State& state, 
					     RealVector& dp)
{
  dp[0] = 1.0;
}
//////////////////////////////////////////////////////////////////////////////

bool Euler2DPuvt::isValid(const RealVector& data)
{
  enum index {P, Ux, Uy, temp};
  const CFreal gamma = getModel()->getGamma();
  const CFreal p = data[P];
  const CFreal T = data[temp];
  const CFreal a = sqrt(gamma*getModel()->getR()*T);
  
  if( ( p < 0.) || (T < 0.) || (a < 0.) ){
    return false;
  }
  
  return true;
}

//////////////////////////////////////////////////////////////////////////////

   } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
