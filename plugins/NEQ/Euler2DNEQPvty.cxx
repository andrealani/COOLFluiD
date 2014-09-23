#include "NEQ.hh"
#include "Euler2DNEQPvty.hh"
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

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQPvty, ConvectiveVarSet, NEQModule, 1>
euler2DNEQPvtyProvider("Euler2DNEQPvty");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPvty::Euler2DNEQPvty(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term),
  _library(CFNULL),
  _dhe(3),
  _ye()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  vector<std::string> names(4 + nbSpecies);
  names[0] = "p";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";

  // Names for the species
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[4 + ie] = "y" + StringOps::to_str(ie);
  }
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPvty::~Euler2DNEQPvty()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPvty::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQPvty::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPvty::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::computePhysicalData(const State& state,
					 RealVector& data)
{
  CFreal p = state[0];
  CFreal T = state[3];

  const RealVector& refData = getModel()->getReferencePhysicalData();

  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*getModel()->getTempRef();

  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // Set the species
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = state[4 + ie];
  }

  _library->setSpeciesFractions(_ye);

  const EquationSubSysDescriptor& eqSS =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint iEqSS = eqSS.getEqSS();

  data[EulerTerm::P] = p;
  data[EulerTerm::T] = T;

  if (iEqSS == 0) {
    _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);

    CFreal rhoDim = _dhe[0];
    data[EulerTerm::RHO] = rhoDim/refData[EulerTerm::RHO];

    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal V2 = u*u + v*v;
    data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
    data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
    _library->frozenGammaAndSoundSpeed(Tdim, pdim, rhoDim, 
				       data[EulerTerm::GAMMA], 
				       data[EulerTerm::A], CFNULL);
    data[EulerTerm::V] = sqrt(V2);
    data[EulerTerm::VX] = u;
    data[EulerTerm::VY] = v;
  }

  if (iEqSS == 1) {
    // here we reuse the previous composition: we neglect the
    // effect of Ye perturbation on the density rho
    data[EulerTerm::RHO] = _library->density(Tdim,pdim)/refData[EulerTerm::RHO];

    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal V2 = u*u + v*v;
    data[EulerTerm::V] = sqrt(V2);
    data[EulerTerm::VX] = state[1];
    data[EulerTerm::VY] = state[2];
  }

  const CFuint nbEqs = eqSS.getNbEqsSS();
  if ((nbEqs == 4 + nbSpecies) || (iEqSS == 1)) {
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      data[firstSpecies + ie] = state[4 + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::computeStateFromPhysicalData(const RealVector& data,
						 State& state)
{
  state[0] = data[EulerTerm::P];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::T];

  // Set the species
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    state[4 + ie] = data[firstSpecies + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DNEQPvty::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::P];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*getModel()->getTempRef();

  // Set the species
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    result[4 + ie] = state[4 + ie];
  }

}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::P];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/(getModel()->getTempRef());

  // Set the species
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    result[4 + ie] = state[4 + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  const CFreal u = state[1]*refData[EulerTerm::V];
  const CFreal v = state[2]*refData[EulerTerm::V];
  CFreal p = state[0]*refData[EulerTerm::P];
  CFreal T = state[3]*getModel()->getTempRef();

  result[0] = p;
  result[1] = u;
  result[2] = v;
  result[3] = T;

  // Set the species
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = result[4 + ie] = state[4 + ie];
  }

  _library->setSpeciesFractions(_ye);

  // density, total enthalpy, Mach, Xi*nbSpecies
  extra.resize(3);
  
  const CFreal V2 = u*u + v*v;
  _library->setDensityEnthalpyEnergy(T, p,_dhe);
  
  extra[0] = _dhe[0];
  extra[1] = _dhe[1] + 0.5*V2;
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  _library->frozenGammaAndSoundSpeed(T,p,extra[0], gamma, a, CFNULL);
  extra[2] = sqrt(V2)/a;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DNEQPvty::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());

  vector<std::string> names(3);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::setup()
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

void Euler2DNEQPvty::computePerturbedPhysicalData
(const State& state, const RealVector& bData, 
 RealVector& data, CFuint iVar) 
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());

  if (iVar == 0 || iVar >= 3) { // p + dp // T + dT // Y + dY
    // this is the case in which the species composition X(p,t,Y) varies
    computePhysicalData(state, data);
  }
  else if(iVar == 1 || iVar == 2) { // u + du // v + dv
    const CFreal u = state[1];
    const CFreal v = state[2];
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
    
    // Terms for species
    const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
    const CFuint nbSpecies = getModel()->getNbScalarVars(0);
    for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
      data[firstSpecies + ie]= bData[firstSpecies + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvty::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  velIDs.resize(2); 
  velIDs[XX] = 1; velIDs[YY] = 2;   
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
