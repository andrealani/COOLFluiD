#include "LTE.hh"
#include "Euler2DLTEDemixPuvt.hh"
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

Environment::ObjectProvider<Euler2DLTEDemixPuvt, ConvectiveVarSet, LTEModule, 1>
euler2DLTEDemixPuvtProvider("Euler2DLTEDemixPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DLTEDemixPuvt::Euler2DLTEDemixPuvt(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term),
  _library(CFNULL),
  _dhe(3),
  _x(),
  _ye()
{
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  
  vector<std::string> names(4 + nbElements);
  names[0] = (!getModel()->isIncompressible()) ? "p" : "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  
  // Names for the elemental fractions
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    names[4 + ie] = "ye" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DLTEDemixPuvt::~Euler2DLTEDemixPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixPuvt::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixPuvt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DLTEDemixPuvt::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixPuvt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::computePhysicalData(const State& state,
					                          RealVector& data)
{  
  const RealVector& refData = getModel()->getReferencePhysicalData();
  
  CFreal T = state[3];
  CFreal pdim = getModel()->getPressureFromState(state[0])*refData[EulerTerm::P];
  CFreal Tdim = T*getModel()->getTempRef();
  
  const CFuint nbEqs =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();

  const CFuint firstElement = getModel()->getFirstScalarVar(0);
  const CFuint nbElements = getModel()->getNbScalarVars(0);

  // Set the elemental fractions
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    _ye[ie] = state[4 + ie];
  }

  _library->setElemFractions(_ye);

  // set the composition
  _library->setComposition(Tdim,pdim);

  data[EulerTerm::P] = state[0];
  data[EulerTerm::T] = T;

  if (nbEqs >= 4) {
    _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);

    CFreal rhoDim = _dhe[0];
    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal V2 = u*u + v*v;

    CFreal gamma = 0.0;
    CFreal adim = 0.0;
    _library->gammaAndSoundSpeed(Tdim,pdim,rhoDim,gamma,adim);

    data[EulerTerm::RHO] = rhoDim/refData[EulerTerm::RHO];
    data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
    data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
    data[EulerTerm::A] = adim/refData[EulerTerm::A];
    data[EulerTerm::V] = sqrt(V2);
    data[EulerTerm::VX] = u;
    data[EulerTerm::VY] = v;
    data[EulerTerm::GAMMA] = gamma;
  }

  if (nbEqs < 4) {
    // here we reuse the previous composition: we neglect the
    // effect of Ye perturbation on the density rho
    data[EulerTerm::RHO] = _library->density(Tdim,pdim)/refData[EulerTerm::RHO];

    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal V2 = u*u + v*v;
    data[EulerTerm::VX] = state[1];
    data[EulerTerm::VY] = state[2];
    data[EulerTerm::V] = sqrt(V2);
  }

  if (nbEqs != 4) {
    for (CFuint ie = 0; ie < nbElements; ++ie) {
      data[firstElement + ie] = state[4 + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::computeStateFromPhysicalData(const RealVector& data,
						  State& state)
{
  state[0] = data[EulerTerm::P];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::T];

  // Set the elemental fractions
  const CFuint firstElement = getModel()->getFirstScalarVar(0);
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    state[4 + ie] = data[firstElement + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DLTEDemixPuvt::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::P];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*getModel()->getTempRef();

  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[4 + ie] = state[4 + ie];
  }

}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::P];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/(getModel()->getTempRef());

  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[4 + ie] = state[4 + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::setDimensionalValuesPlusExtraValues
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

vector<std::string> Euler2DLTEDemixPuvt::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());
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

void Euler2DLTEDemixPuvt::setup()
{
  MultiScalarVarSet<Euler2DVarSet>::setup();

  // set the equation set data for each of the equation subsets
  // first equation subset
  Euler2DVarSet::getEqSetData().resize(1);
  Euler2DVarSet::getEqSetData()[0].setup(0,0,4);
  
  // second equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData().resize(1);
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].setup
    (1,4,getModel()->getNbScalarVars(0));
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());

  _ye.resize(getModel()->getNbScalarVars(0));
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvt::computePerturbedPhysicalData(const Framework::State& currState,
						    const RealVector& bData,
						    RealVector& data,
						    CFuint iVar)
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());
  
  if (iVar == 0 || iVar >= 3) { // p + dp // T + dT // Y + dY
    computePhysicalData(currState, data);
  }
  else if(iVar == 1 || iVar == 2) { // u + du // v + dv
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
    data[EulerTerm::GAMMA] = bData[EulerTerm::GAMMA];
    data[EulerTerm::E] = bData[EulerTerm::E] +
      0.5* (V2 - bData[EulerTerm::V]*bData[EulerTerm::V]);
    
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
