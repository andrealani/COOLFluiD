#include "NavierStokes/NavierStokes.hh"
#include "Euler3DRoe.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DRoe, ConvectiveVarSet, NavierStokesModule, 1> 
euler3DRoeProvider("Euler3DRoe");

//////////////////////////////////////////////////////////////////////////////

Euler3DRoe::Euler3DRoe(Common::SafePtr<BaseTerm> term) :
  Euler3DVarSet(term)
{
  vector<std::string> names(5);
  names[0] = "z0";
  names[1] = "z1";
  names[2] = "z2";
  names[3] = "z3";
  names[4] = "z4";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler3DRoe::~Euler3DRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRoe::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler3DRoe::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRoe::computeEigenValuesVectors(RealMatrix& rightEv,
          RealMatrix& leftEv,
          RealVector& eValues,
          const RealVector& normal)
{
 throw Common::NotImplementedException (FromHere(),"Euler3DRoe::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DRoe::getBlockSeparator() const
{
  return 5;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRoe::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DRoe::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRoe::computePhysicalData(const State& state, RealVector& data)
  
{
  const CFreal sqRho = state[0];
  const CFreal rho = sqRho*sqRho;
  const CFreal u = state[1]/sqRho;
  const CFreal v = state[2]/sqRho;
  const CFreal w = state[3]/sqRho;
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal gamma = getModel()->getGamma();
  const CFreal p = (gamma - 1.)/gamma*sqRho*(state[4] - 0.5*V2*sqRho);
  const CFreal pOvRho = p/rho;
  
  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = p;
  data[EulerTerm::H] = state[4]/sqRho;
  data[EulerTerm::E] = data[EulerTerm::H] - pOvRho;
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = pOvRho/getModel()->getR();
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::VZ] = w;
  data[EulerTerm::GAMMA] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRoe::computeStateFromPhysicalData(const RealVector& data,
				     State& state)
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);
  state[0] = sqRho;
  state[1] = sqRho*data[EulerTerm::VX];
  state[2] = sqRho*data[EulerTerm::VY];
  state[3] = sqRho*data[EulerTerm::VZ];
  state[4] = sqRho*data[EulerTerm::H];
}
      
//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DRoe::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2] + state[3]*state[3])/state[0];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler3DRoe::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DRoe::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRoe::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DRoe::setAdimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
