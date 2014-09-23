#include "LTE.hh"
#include "Euler3DLTEDemixPvt.hh"
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

Environment::ObjectProvider<Euler3DLTEDemixPvt, ConvectiveVarSet, LTEModule, 1>
euler3DLTEDemixPvtProvider("Euler3DLTEDemixPvt");

//////////////////////////////////////////////////////////////////////////////

Euler3DLTEDemixPvt::Euler3DLTEDemixPvt(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler3DVarSet>(term),
  _library(CFNULL),
  _dhe(3),
  _x(),
  _ye()
{
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  
  vector<std::string> names(5 + nbElements);
  names[0] = (!getModel()->isIncompressible()) ? "p" : "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "T";
  
  // Names for the elemental fractions
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    names[5 + ie] = "ye" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler3DLTEDemixPvt::~Euler3DLTEDemixPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixPvt::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::computeEigenValuesVectors(RealMatrix& rightEv,
					   RealMatrix& leftEv,
					   RealVector& eValues,
					   const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixPvt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DLTEDemixPvt::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::splitJacobian(RealMatrix& jacobPlus,
				RealMatrix& jacobMin,
				RealVector& eValues,
				const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixPvt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::computePhysicalData(const State& state,
					     RealVector& data)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  
  CFreal pdim = getModel()->getPressureFromState(state[0])*refData[EulerTerm::P];
  CFreal T = state[4];
  CFreal Tdim = T*getModel()->getTempRef();
  
  // Set the elemental fractions
  const CFuint firstElement = getModel()->getFirstScalarVar(0);
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    _ye[ie] = data[firstElement + ie] = state[5 + ie];
  }
  
  _library->setElemFractions(_ye);

  // set the composition
  _library->setComposition(Tdim,pdim);
  _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);

  CFreal rhoDim = _dhe[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal V2 = u*u + v*v +w*w;
  
  CFreal gamma = 0.0;
  CFreal adim = 0.0;
  _library->gammaAndSoundSpeed(Tdim,pdim,rhoDim,gamma,adim);
  
  data[EulerTerm::RHO] = rhoDim/refData[EulerTerm::RHO];
  data[EulerTerm::P] = state[0];
  data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
  data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
  data[EulerTerm::A] = adim/refData[EulerTerm::A];
  data[EulerTerm::T] = T;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::VZ] = w;
  data[EulerTerm::GAMMA] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::computeStateFromPhysicalData(const RealVector& data,
						 State& state)
{
  state[0] = data[EulerTerm::P];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::VZ];
  state[4] = data[EulerTerm::T];

  // Set the elemental fractions
  const CFuint firstElement = getModel()->getFirstScalarVar(0);
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    state[5 + ie] = data[firstElement + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DLTEDemixPvt::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2] + state[3]*state[3]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::setDimensionalValues(const State& state,
					      RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  
  result[0] = state[0]*refData[EulerTerm::P];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::V];
  result[4] = state[4]*getModel()->getTempRef();
  
  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[5 + ie] = state[5 + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::setAdimensionalValues(const State& state,
					  RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  
  result[0] = state[0]/refData[EulerTerm::P];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/refData[EulerTerm::V];
  result[4] = state[4]/(getModel()->getTempRef());
  
  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[5 + ie] = state[5 + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  const CFreal u = state[1]*refData[EulerTerm::V];
  const CFreal v = state[2]*refData[EulerTerm::V];
  const CFreal w = state[3]*refData[EulerTerm::V];
  CFreal p = getModel()->getPressureFromState(state[0])*refData[EulerTerm::P];
  CFreal T = state[4]*getModel()->getTempRef();
  
  result[0] = state[0];
  result[1] = u;
  result[2] = v;
  result[3] = w;
  result[4] = T;
  
  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    _ye[ie] = result[5 + ie] = state[5 + ie];
  }

  _library->setElemFractions(_ye);

  const CFuint nbSpecies = _library->getNbSpecies();
  // density, total enthalpy, Mach, Xi*nbSpecies
  extra.resize(3 + nbSpecies);
  
  _x.resize(nbSpecies);
  
  // set the composition
  _library->setComposition(T,p,&_x);
  _library->setDensityEnthalpyEnergy(T,p,_dhe);
  
  const CFreal V2 = u*u + v*v + w*w;
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

vector<std::string> Euler3DLTEDemixPvt::getExtraVarNames() const
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

void Euler3DLTEDemixPvt::setup()
{
  MultiScalarVarSet<Euler3DVarSet>::setup();
  
  // set the equation set data for each of the equation subsets
  // first equation subset
  Euler3DVarSet::getEqSetData().resize(1);
  Euler3DVarSet::getEqSetData()[0].setup(0,0,5);
  
  // second equation subset
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData().resize(1);
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData()[0].setup
    (1,5,getModel()->getNbScalarVars(0));
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());

  _ye.resize(getModel()->getNbScalarVars(0));
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvt::computePerturbedPhysicalData(const Framework::State& currState,
						  const RealVector& bData,
						  RealVector& data,
						  CFuint iVar)
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());

  if (iVar == 0 || iVar >= 4) { // p + dp // T + dT // Y + dY
    computePhysicalData(currState, data);
  }
  else if(iVar == 1 || iVar == 2 || iVar == 3) { // u + du // v + dv // w + dw
    const CFreal u = currState[1];
    const CFreal v = currState[2];
    const CFreal w = currState[3];
    const CFreal V2 = u*u + v*v + w*w;
    
    data[EulerTerm::VX] = u;
    data[EulerTerm::VY] = v;
    data[EulerTerm::VZ] = w;
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
