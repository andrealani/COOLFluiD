#include "NavierStokes/NavierStokes.hh"
#include <numeric>

#include "Euler1DPvt.hh"
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

Environment::ObjectProvider<Euler1DPvt, ConvectiveVarSet, NavierStokesModule, 1> euler1DPvtProvider("Euler1DPvt");

//////////////////////////////////////////////////////////////////////////////

Euler1DPvt::Euler1DPvt(Common::SafePtr<BaseTerm> term) :
  Euler1DVarSet(term)
{
  vector<std::string> names(3);
  names[0] = (!getModel()->isIncompressible()) ? "p" : "dp";
  names[1] = "v";
  names[2] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler1DPvt::~Euler1DPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler1DPvt::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::computeEigenValuesVectors(RealMatrix& rightEv,
                                        RealMatrix& leftEv,
                                        RealVector& eValues,
                                        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler1DPvt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler1DPvt::getBlockSeparator() const
{
  return 3;
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler1DPvt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::computePhysicalData(const State& state, RealVector& data)
{
  const CFreal p = getModel()->getPressureFromState(state[0]);
  cf_assert(p > 0);
  const CFreal u = state[1];
  const CFreal T = state[2];
  const CFreal rho = getModel()->getDensity(p,T);
  const CFreal V2 = u*u;
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
  data[EulerTerm::GAMMA] = gamma;
  
  if (this->_physDataNeedCoordinates) {
    const RealVector& node = state.getCoordinates();
    data[EulerTerm::XP] = node[XX];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::computeStateFromPhysicalData(const RealVector& data, State& state)
{
  state[0] = data[EulerTerm::P]; // p or dp
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler1DPvt::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::P];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::P];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/(getModel()->getTempRef());
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DPvt::computePressureDerivatives(const Framework::State& state, 
					     RealVector& dp)
{
  dp[0] = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

bool Euler1DPvt::isValid(const RealVector& data)
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
