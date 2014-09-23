#include "NavierStokes/NavierStokes.hh"
#include <numeric>

#include "Euler3DPrim.hh"
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

Environment::ObjectProvider<Euler3DPrim, ConvectiveVarSet, NavierStokesModule, 1> Euler3DPrimProvider("Euler3DPrim");

//////////////////////////////////////////////////////////////////////////////

Euler3DPrim::Euler3DPrim(Common::SafePtr<BaseTerm> term) :
  Euler3DVarSet(term)
{
   vector<std::string> names(5);
   names[0] = "rho";
   names[1] = "u";
   names[2] = "v";
   names[3] = "w";
   names[4] = (!getModel()->isIncompressible()) ? "p" : "dp";
   setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler3DPrim::~Euler3DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrim::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler3DPrim::computeJacobians()");

}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrim::computeEigenValuesVectors(RealMatrix& rightEv,
          RealMatrix& leftEv,
          RealVector& eValues,
          const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DPrim::computeEigenValuesVectors()");

}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DPrim::getBlockSeparator() const
{
  return 5;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrim::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DPrim::splitJacobian()");

}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrim::computePhysicalData(const State& state, RealVector& data)
  
{
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal p = getModel()->getPressureFromState(state[4]);
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal pOvRho = p/rho;

  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = state[4];
  data[EulerTerm::H] = gammaDivGammaMinus1*pOvRho + 0.5*V2;
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

void Euler3DPrim::computeStateFromPhysicalData(const RealVector& data,
					  State& state)
{
  state[0] = data[EulerTerm::RHO];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::VZ];
  state[4] = data[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DPrim::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2] + state[3]*state[3]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrim::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::RHO];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::V];
  result[4] = state[4]*refData[EulerTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPrim::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::RHO];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/refData[EulerTerm::V];
  result[4] = state[4]/refData[EulerTerm::P];
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
