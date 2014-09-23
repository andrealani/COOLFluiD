#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRoe.hh"
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

Environment::ObjectProvider<Euler2DRoe, ConvectiveVarSet, NavierStokesModule, 1> 
euler2DRoeProvider("Euler2DRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DRoe::Euler2DRoe(Common::SafePtr<BaseTerm> term) :
  Euler2DVarSet(term)
{
  vector<std::string> names(4);
  names[0] = "z0";
  names[1] = "z1";
  names[2] = "z2";
  names[3] = "z3";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRoe::~Euler2DRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoe::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler2DRoe::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoe::computeEigenValuesVectors(RealMatrix& rightEv,
          RealMatrix& leftEv,
          RealVector& eValues,
          const RealVector& normal)
{
 throw Common::NotImplementedException (FromHere(),"Euler2DRoe::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DRoe::getBlockSeparator() const
{
  return 4;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoe::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DRoe::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoe::computePhysicalData(const State& state, RealVector& data)
  
{
  const CFreal sqRho = state[0];
  const CFreal rho = sqRho*sqRho;
  const CFreal u = state[1]/sqRho;
  const CFreal v = state[2]/sqRho;
  const CFreal V2 = u*u + v*v;
  const CFreal gamma = getModel()->getGamma();
  const CFreal p = (gamma - 1.)/gamma*sqRho*(state[3] - 0.5*V2*sqRho);
  const CFreal pOvRho = p/rho;
  
  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = p;
  data[EulerTerm::H] = state[3]/sqRho;
  data[EulerTerm::E] = data[EulerTerm::H] - pOvRho;
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = pOvRho/getModel()->getR();
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::GAMMA] = gamma;
  
  //  cout << "data  = " << data << endl;
  //cout << "state = " << state << endl << endl;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoe::computeStateFromPhysicalData(const RealVector& data,
					      State& state)
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);
  state[0] = sqRho;
  state[1] = sqRho*data[EulerTerm::VX];
  state[2] = sqRho*data[EulerTerm::VY];
  state[3] = sqRho*data[EulerTerm::H];
}
      
      //////////////////////////////////////////////////////////////////////////////

CFreal Euler2DRoe::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2])/state[0];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DRoe::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DRoe::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoe::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DRoe::setAdimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
