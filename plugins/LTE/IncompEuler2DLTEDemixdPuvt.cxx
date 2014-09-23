#include "LTE.hh"
#include "IncompEuler2DLTEDemixdPuvt.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

/////////////////////////////////////////////////////////////////////
/////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<IncompEuler2DLTEDemixdPuvt, ConvectiveVarSet, LTEModule, 1>
incompEuler2DLTEDemixdPuvtProvider("IncompEuler2DLTEDemixdPuvt");

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DLTEDemixdPuvt::IncompEuler2DLTEDemixdPuvt(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<IncompEuler2DVarSet>(term),
  _library(CFNULL),
  _dhe(3),
  _x(),
  _ye()
{
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  
  vector<std::string> names(4 + nbElements);
  names[0] = "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";

 // Names for the elemental fractions
 for (CFuint ie = 0; ie < nbElements; ++ie) {
   names[4 + ie] = "ye" + StringOps::to_str(ie);
 }

 setVarNames(names);
}

/////////////////////////////////////////////////////////////////////
/////////

IncompEuler2DLTEDemixdPuvt::~IncompEuler2DLTEDemixdPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

CFuint IncompEuler2DLTEDemixdPuvt::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixdPuvt::computePhysicalData(const State& state, RealVector& data)

{
  CFreal dp = state[0];
  CFreal T  = state[3];
  CFreal p  = getModel()->getThermodynamPressInf() + dp;
  
  const RealVector& refData = getModel()->getReferencePhysicalData();

  CFreal pdim = p; 
  CFreal Tdim = T*refData[IncompEulerTerm::T];

  const CFuint nbEqs = PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();
  const CFuint firstElement = getModel()->getFirstScalarVar(0);  
  const CFuint nbElements = getModel()->getNbScalarVars(0);

  // Set the elemental fractions
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    _ye[ie] = state[4 + ie];
  }
  
  _library->setElemFractions(_ye);
  

  // set the composition
  
  _library->setComposition(Tdim,pdim);

  data[IncompEulerTerm::dP]  = dp;
  data[IncompEulerTerm::T] = T;

  if (nbEqs >= 4) {  
    _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);

    CFreal rhoDim = _dhe[0];
    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal V2 = u*u + v*v; 

    CFreal gamma = 0.0;
    CFreal adim = 0.0;
    _library->gammaAndSoundSpeed(Tdim,pdim,rhoDim,gamma,adim);

    data[IncompEulerTerm::RHO] = rhoDim/refData[IncompEulerTerm::RHO];
    data[IncompEulerTerm::H] = _dhe[1]/refData[IncompEulerTerm::H] + 0.5*V2;
    data[IncompEulerTerm::E] = _dhe[2]/refData[IncompEulerTerm::H] + 0.5*V2;
    data[IncompEulerTerm::A] = adim/refData[IncompEulerTerm::A];
    data[IncompEulerTerm::V] = sqrt(V2);
    data[IncompEulerTerm::VX] = u;
    data[IncompEulerTerm::VY] = v;
  }

  if (nbEqs < 4) {
    // here we reuse the previous composition: we neglect the
    // effect of Ye perturbation on the density rho
    
    data[IncompEulerTerm::RHO] = _library->density(Tdim,pdim)/refData[IncompEulerTerm::RHO];
    
    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal V2 = u*u + v*v;
    data[IncompEulerTerm::VX] = state[1];
    data[IncompEulerTerm::VY] = state[2];
    data[IncompEulerTerm::V] = sqrt(V2);
  }
  
  if (nbEqs != 4) {
    for (CFuint ie = 0; ie < nbElements; ++ie) {
      data[firstElement + ie] = state[4 + ie];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixdPuvt::computeStateFromPhysicalData(const RealVector& data, State& state)
{
  state[0] = data[IncompEulerTerm::dP];
  state[1] = data[IncompEulerTerm::VX];
  state[2] = data[IncompEulerTerm::VY];
  state[3] = data[IncompEulerTerm::T];

  // Set the elemental fractions
  const CFuint firstElement = getModel()->getFirstScalarVar(0);
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    state[4 + ie] = data[firstElement + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal IncompEuler2DLTEDemixdPuvt::getSpeed(const State& state) const
{
   return sqrt(state[1]*state[1] + state[2]*state[2]);
}


//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixdPuvt::setDimensionalValues
(const State& state,RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[IncompEulerTerm::dP];
  result[1] = state[1]*refData[IncompEulerTerm::V];
  result[2] = state[2]*refData[IncompEulerTerm::V];
  result[3] = state[3]*getModel()->getTempRef();

  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[4 + ie] = state[4 + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixdPuvt::setAdimensionalValues
(const Framework::State& state,RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[IncompEulerTerm::dP];
  result[1] = state[1]/refData[IncompEulerTerm::V];
  result[2] = state[2]/refData[IncompEulerTerm::V];
  result[3] = state[3]/(getModel()->getTempRef());

  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[4 + ie] = state[4 + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixdPuvt::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result, RealVector& extra)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal u = state[1]*refData[IncompEulerTerm::V];
  const CFreal v = state[2]*refData[IncompEulerTerm::V];
  const CFreal dp = state[0]*refData[IncompEulerTerm::dP];
  CFreal T = state[3]*getModel()->getTempRef();
  
  result[0] = dp;
  result[1] = u;
  result[2] = v;
  result[3] = T;

  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    _ye[ie] = result[4 + ie] = state[4 + ie];
  }

  _library->setElemFractions(_ye);
  const CFuint nbSpecies = _library->getNbSpecies();

  // density, total enthalpy, Mach, Xi*nbSpecies
  extra.resize(3 + nbSpecies);
  
  _x.resize(nbSpecies);
  
  // TODO This is not correct in adimensional case ???
  CFreal p = getModel()->getThermodynamPressInf() + dp;

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

vector<std::string> IncompEuler2DLTEDemixdPuvt::getExtraVarNames() const
{
  cf_assert (_library.isNotNull());
  const CFuint nbSpecies = _library->getNbSpecies();

  vector<std::string> names(3 + nbSpecies);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  for (CFuint is = 0; is < nbSpecies; ++is) {
    names[3 + is] = "xc" + StringOps::to_str(is);
  }

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixdPuvt::setup()
{
  MultiScalarVarSet<IncompEuler2DVarSet>::setup();

  // set the equation set data for each of the equation subsets
  // first equation subset
  IncompEuler2DVarSet::getEqSetData().resize(1);
  IncompEuler2DVarSet::getEqSetData()[0].setup(0,0,4);

  // second equation subset
  MultiScalarVarSet<IncompEuler2DVarSet>::getEqSetData().resize(1);
  MultiScalarVarSet<IncompEuler2DVarSet>::getEqSetData()[0].setup(1,4,getModel()->getNbScalarVars(0));

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert (_library.isNotNull());

  _ye.resize(getModel()->getNbScalarVars(0));
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixdPuvt::computePerturbedPhysicalData(const Framework::State& currState,
						    const RealVector& bData,
						    RealVector& data,
						    CFuint iVar)
{
  cf_assert (iVar < PhysicalModelStack::getActive()->getNbEq());
  if (iVar == 0 || iVar >= 3) { // p + dp // T + dT // Y + dY
    computePhysicalData(currState, data);
  }
  else if(iVar == 1 || iVar == 2) { // u + du // v + dv
    const CFreal u = currState[1];
    const CFreal v = currState[2];
    const CFreal V2 = u*u + v*v;
    
    data[IncompEulerTerm::VX] = u;
    data[IncompEulerTerm::VY] = v;
    data[IncompEulerTerm::RHO] = bData[IncompEulerTerm::RHO];
    data[IncompEulerTerm::V] = sqrt(V2);
    data[IncompEulerTerm::dP] = bData[IncompEulerTerm::dP];
    data[IncompEulerTerm::H] = bData[IncompEulerTerm::H] + 
      0.5* (V2 - bData[IncompEulerTerm::V]*bData[IncompEulerTerm::V]);
    
    data[IncompEulerTerm::A] = bData[IncompEulerTerm::A];
    data[IncompEulerTerm::E] = bData[IncompEulerTerm::E] +
      0.5* (V2 - bData[IncompEulerTerm::V]*bData[IncompEulerTerm::V]);
    
    // Terms for elemental fractions
    const CFuint firstElement = getModel()->getFirstScalarVar(0);
    const CFuint nbElements = getModel()->getNbScalarVars(0);
    for (CFuint ie = 0;  ie < nbElements; ++ie) {
      data[firstElement + ie]= bData[firstElement + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




