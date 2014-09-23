#include "LTE.hh"
#include "Euler2DPuvtLTE.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DPuvtLTE, ConvectiveVarSet, LTEModule, 1>
euler2DPuvtLTEProvider("Euler2DPuvtLTE");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtLTE::Euler2DPuvtLTE(Common::SafePtr<BaseTerm> term) :
  Euler2DVarSet(term),
  _library(CFNULL),
  _dhe(3),
  _x()
{
  vector<std::string> names(4);
  names[0] = (!getModel()->isIncompressible()) ? "p" : "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtLTE::~Euler2DPuvtLTE()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler2DPuvtLTE::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DPuvtLTE::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DPuvtLTE::getBlockSeparator() const
{
  return 4;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DPuvtLTE::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::computePhysicalData(const State& state,
					 RealVector& data)
{ 
  const RealVector& refData = getModel()->getReferencePhysicalData();
    
  CFreal pdim = getModel()->getPressureFromState(state[0])*refData[EulerTerm::P];
  // correction in case pdim=0 (e.g. when starting from scratch in compressible mode)
  pdim = (pdim > 0.01) ? pdim :getModel()->getPressInfComp()*refData[EulerTerm::P];
  cf_assert(pdim > 0.);
    
  CFreal T = state[3];
  CFreal Tdim = T*getModel()->getTempRef();
  cf_assert(_library.isNotNull());
 
  CFLog(DEBUG_MAX, "Euler2DPuvtLTE::computePhysicalData() => Tdim = " << Tdim << "\n");
  CFLog(DEBUG_MAX, "Euler2DPuvtLTE::computePhysicalData() => pdim = " << pdim << "\n");   
  // set the composition
  _library->setComposition(Tdim,pdim); 
  _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
  CFLog(DEBUG_MAX, "Euler2DPuvtLTE::computePhysicalData() => _dhe= " << _dhe  << "\n");

  CFreal rhoDim = _dhe[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal V2 = u*u + v*v;

  CFreal gamma = 0.0;
  CFreal adim = 0.0;
  _library->gammaAndSoundSpeed(Tdim,pdim,rhoDim,gamma,adim);
  CFLog(DEBUG_MAX, "Euler2DPuvtLTE::computePhysicalData() => gamma = " << gamma <<  "\n");
  CFLog(DEBUG_MAX, "Euler2DPuvtLTE::computePhysicalData() => adim = " << adim <<  "\n");

  data[EulerTerm::RHO] = rhoDim/refData[EulerTerm::RHO];
  data[EulerTerm::P] = state[0];
  data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
  data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
  data[EulerTerm::A] = adim/refData[EulerTerm::A];
  data[EulerTerm::T] = T;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::T] = T;
  data[EulerTerm::GAMMA] = gamma;
  
  CFLog(DEBUG_MAX, "Euler2DPuvtLTE::computePhysicalData() => data = " << data << "\n");
  CFLog(DEBUG_MAX, "Euler2DPuvtLTE::computePhysicalData() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::computeStateFromPhysicalData
(const RealVector& data, State& state)
{
  state[0] = data[EulerTerm::P];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DPuvtLTE::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::P];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::setAdimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::P];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/(getModel()->getTempRef());
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  const CFreal u = state[1]*refData[EulerTerm::V];
  const CFreal v = state[2]*refData[EulerTerm::V];
  CFreal p = getModel()->getPressureFromState(state[0])*refData[EulerTerm::P];
  CFreal T = state[3]*getModel()->getTempRef();

  result[0] = state[0];
  result[1] = u;
  result[2] = v;
  result[3] = T;

  const CFuint nbSpecies = _library->getNbSpecies();
  // density, total enthalpy, Mach, Xi*nbSpecies
  extra.resize(3 + nbSpecies);

  _x.resize(_library->getNbSpecies());

  // set the composition
  _library->setComposition(T,p,&_x);
  const CFreal V2 = u*u + v*v;
  _library->setDensityEnthalpyEnergy(T,p,_dhe);
  extra[0] = _dhe[0];
  extra[1] = _dhe[1] + 0.5*V2;
  
  CFLog(VERBOSE, "Euler2DPuvtLTE::setDimensionalValuesPlusExtraValues() => dhe = " << _dhe << "\n");
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  _library->gammaAndSoundSpeed(T,p,_dhe[0],gamma,a);
  extra[2] = sqrt(V2)/a;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    extra[3 + is] = _x[is];
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DPuvtLTE::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());
  const CFuint nbSpecies = _library->getNbSpecies();

  vector<std::string> names(3 + nbSpecies);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  for (CFuint is = 0; is < nbSpecies; ++is) {
    names[3 + is] = "xc" + StringOps::to_str(is);;
  }

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::setup()
{
  //  fluct split has to call setup()
  Euler2DVarSet::setup();

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTE::computePerturbedPhysicalData(const Framework::State& currState,
						  const RealVector& bData,
						  RealVector& data,
						  CFuint iVar)
{
  cf_assert(iVar < 4);

  if (iVar == 0 || iVar == 3) { // p + dp // T + dT
    computePhysicalData(currState, data);
  }
  else { // u + du // v + dv
    const CFreal u = currState[1];
    const CFreal v = currState[2];
    const CFreal V2 = u*u + v*v;

    data[EulerTerm::VX] = u;
    data[EulerTerm::VY] = v;
    data[EulerTerm::RHO] = bData[EulerTerm::RHO];
    data[EulerTerm::V] = sqrt(V2);
    data[EulerTerm::P] = bData[EulerTerm::P];
    data[EulerTerm::H] = bData[EulerTerm::H] +
      0.5* (V2 - bData[EulerTerm::V]*bData[EulerTerm::V]);
    data[EulerTerm::A] = bData[EulerTerm::A];
    data[EulerTerm::E] = bData[EulerTerm::E] +
      0.5* (V2 - bData[EulerTerm::V]*bData[EulerTerm::V]);
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
