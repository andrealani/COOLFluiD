#include "LTE.hh"
#include "IncompEuler2DdPuvtLTE.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////
namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<IncompEuler2DdPuvtLTE, ConvectiveVarSet, LTEModule, 1>
incompEuler2DdPuvtLTEProvider("IncompEuler2DdPuvtLTE");

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtLTE::IncompEuler2DdPuvtLTE(Common::SafePtr<Framework::BaseTerm> term) :
  IncompEuler2DVarSet(term),
  _library(CFNULL),
  _dhe(3),
  _x()
{
  vector<std::string> names(4);
  names[0] = "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtLTE::~IncompEuler2DdPuvtLTE()
{
}

//////////////////////////////////////////////////////////////////////////////

CFuint IncompEuler2DdPuvtLTE::getBlockSeparator() const
{
  return 4;
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTE::computePhysicalData(const State& state,
						RealVector& data)
{
  const CFreal dp = state[0];
  const CFreal u  = state[1];
  const CFreal v  = state[2];
  const CFreal T  = state[3];
  const CFreal p  = getModel()->getThermodynamPressInf() + dp;
  const CFreal V2  = u*u + v*v;

  const RealVector& refData = getModel()->getReferencePhysicalData();

  CFreal pdim = p; // *getModel()->getPressRef(); watch out with this reference pressure
  CFreal Tdim = T*refData[IncompEulerTerm::T];
  _library->setComposition(Tdim,pdim);
  _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
  CFreal rhoDim = _dhe[0];
  CFreal gamma = 0.0;
  CFreal adim = 0.0;
  _library->gammaAndSoundSpeed(Tdim, pdim,rhoDim, gamma, adim);

  data[IncompEulerTerm::VX]  = u;
  data[IncompEulerTerm::VY]  = v;
  data[IncompEulerTerm::RHO] = rhoDim/refData[IncompEulerTerm::RHO];
  data[IncompEulerTerm::V]   = sqrt(V2);
  data[IncompEulerTerm::dP]  = dp;
  data[IncompEulerTerm::T]   = T;
  data[IncompEulerTerm::H]   = _dhe[1]/refData[IncompEulerTerm::H] + 0.5*V2;
  data[IncompEulerTerm::E]   = _dhe[2]/refData[IncompEulerTerm::E] + 0.5*V2;
  data[IncompEulerTerm::A]   = adim/refData[IncompEulerTerm::A];
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTE::computeStateFromPhysicalData(const RealVector& data,
						    State& state)
{
  state[0] = data[IncompEulerTerm::dP];
  state[1] = data[IncompEulerTerm::VX];
  state[2] = data[IncompEulerTerm::VY];
  state[3] = data[IncompEulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal IncompEuler2DdPuvtLTE::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2]);
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTE::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[IncompEulerTerm::dP];
  result[1] = state[1]*refData[IncompEulerTerm::V];
  result[2] = state[2]*refData[IncompEulerTerm::V];
  result[3] = state[3]*refData[IncompEulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTE::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[IncompEulerTerm::dP];
  result[1] = state[1]/refData[IncompEulerTerm::V];
  result[2] = state[2]/refData[IncompEulerTerm::V];
  result[3] = state[3]/refData[IncompEulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> IncompEuler2DdPuvtLTE::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());
  const CFuint nbSpecies = _library->getNbSpecies();

  vector<std::string> names(3 + nbSpecies);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  for (CFuint is = 0; is < nbSpecies; ++is) {
    names[3 + is] = "xc" + Common::StringOps::to_str(is);;
  }

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTE::setup()
{
  IncompEuler2DVarSet::setup();

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTE::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result, RealVector& extra)
{
  setDimensionalValues(state, result);

  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal u = state[1]*refData[IncompEulerTerm::V];
  const CFreal v = state[2]*refData[IncompEulerTerm::V];
  const CFreal dp = state[0]*refData[IncompEulerTerm::dP];
  CFreal T = state[3]*refData[IncompEulerTerm::T];

  const CFuint nbSpecies = _library->getNbSpecies();
  // density, total enthalpy, Mach, Xi*nbSpecies
  extra.resize(3 + nbSpecies);

  _x.resize(_library->getNbSpecies());

  /// @TODO this is not correct in the adimensional case !!!!
  CFreal p = getModel()->getThermodynamPressInf() + dp; // *getModel()->getPressRef();

  // set the composition
  _library->setComposition(T,p,&_x);
  _library->setDensityEnthalpyEnergy(T,p,_dhe);
  const CFreal V2 = u*u + v*v;

  CFreal gamma = 0.0;
  CFreal adim = 0.0;
  _library->gammaAndSoundSpeed(T,p,_dhe[0],gamma,adim);

  extra[0] = _dhe[0];
  extra[1] = _dhe[1] + 0.5*V2;
  extra[2] = sqrt(V2)/adim;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    extra[3 + is] = _x[is];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTE::computePerturbedPhysicalData(const Framework::State& currState,
						  const RealVector& bData,
						  RealVector& data,
						  CFuint iVar)
{
  if (iVar != 1 && iVar != 2) {// p + dp // T + dT // other variables ...
    computePhysicalData(currState, data);
  }
  else if (iVar == 1 || iVar == 2) { // u + du // v + dv
    const CFreal u = currState[1];
    const CFreal v = currState[2];
    const CFreal V2 = u*u + v*v;
    
    data[IncompEulerTerm::VX] = u;
    data[IncompEulerTerm::VY] = v;
    data[IncompEulerTerm::RHO] = bData[IncompEulerTerm::RHO];
    data[IncompEulerTerm::V] = sqrt(V2);
    data[IncompEulerTerm::dP] = bData[IncompEulerTerm::dP];
    data[IncompEulerTerm::T] = bData[IncompEulerTerm::T];
    data[IncompEulerTerm::H] = bData[IncompEulerTerm::H] +
      0.5* (V2 - bData[IncompEulerTerm::V]*bData[IncompEulerTerm::V]);
    data[IncompEulerTerm::E] = bData[IncompEulerTerm::E] +
      0.5* (V2 - bData[IncompEulerTerm::V]*bData[IncompEulerTerm::V]);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
