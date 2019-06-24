#include "NavierStokes/NavierStokes.hh"
#include <numeric>

#include "NavierStokes/Euler3DRhovt.hh"
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

Environment::ObjectProvider<Euler3DRhovt, ConvectiveVarSet, NavierStokesModule, 1>
euler3DRhovtProvider("Euler3DRhovt");

//////////////////////////////////////////////////////////////////////////////

Euler3DRhovt::Euler3DRhovt(Common::SafePtr<BaseTerm> term) :
  Euler3DVarSet(term)
{
  vector<std::string> names(5);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "T";
  setVarNames(names);
}
      
//////////////////////////////////////////////////////////////////////////////

Euler3DRhovt::~Euler3DRhovt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler3DRhovt::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::computeEigenValuesVectors(RealMatrix& rightEv,
                                        RealMatrix& leftEv,
                                        RealVector& eValues,
                                        const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DRhovt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DRhovt::getBlockSeparator() const
{
  return 5;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DRhovt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::computePhysicalData(const State& state, RealVector& data)
{
  // CFLogInfo("state = " << state << "\n");
  
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal T = state[4];
  const CFreal R = getModel()->getR();
  const CFreal p = rho*R*T;
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal pOvRho = p/rho;
  
  data[EulerTerm::RHO] = rho;
  data[EulerTerm::P] = p;
  data[EulerTerm::H] = gammaDivGammaMinus1*pOvRho + 0.5*V2;
  data[EulerTerm::E] = data[EulerTerm::H] - pOvRho;
  data[EulerTerm::A] = sqrt(gamma*pOvRho);
  data[EulerTerm::T] = T;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::VZ] = w;
  data[EulerTerm::GAMMA] = gamma;
  
  // CFLogInfo("data = " << data << "\n\n");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::computeStateFromPhysicalData(const RealVector& data,
				       State& state)
{
  const CFreal R = getModel()->getR();
  state[0] = data[EulerTerm::RHO];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::VZ];
  state[4] = data[EulerTerm::P]/(R*data[EulerTerm::RHO]);
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DRhovt::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2] + state[3]*state[3]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::RHO];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::V];
  result[4] = state[4]*getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  
  result[0] = state[0]/refData[EulerTerm::RHO];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/refData[EulerTerm::V];
  result[4] = state[4]/(getModel()->getTempRef());
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DRhovt::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  // first set the state variables
  Euler3DRhovt::setDimensionalValues(state,result);

  const CFreal rho = result[0];
  const CFreal T   = result[4];
  const CFreal R = getModel()->getR();
  const CFreal p = rho*R*T;
  const CFreal u = result[1];
  const CFreal v = result[2];
  const CFreal w = result[3];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal pOvRho = p/rho;
  const CFreal M = sqrt(V2/(gamma*pOvRho));
  
  extra.resize(4);
  extra[0] = rho;
  extra[1] = gammaDivGammaMinus1*pOvRho + 0.5*V2;
  extra[2] = M;
  extra[3] = p;
  
  /*
    static CFreal minM = 1e10; if (M<minM) minM = M; 
    static CFreal maxM = 0.;   if (M>maxM) maxM = M;
    CFLog(INFO, "Mach [min, max] = [" << minM << ", " << maxM << "]\n");
  */
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler3DRhovt::getExtraVarNames() const
{  
  vector<std::string> names(4);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";
  names[3] = "p";
  
  return names;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
