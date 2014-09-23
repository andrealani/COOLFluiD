#include "NEQ.hh"
#include "Euler3DNEQPvty.hh"
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

Environment::ObjectProvider<Euler3DNEQPvty, ConvectiveVarSet, NEQModule, 1>
euler3DNEQPvtyProvider("Euler3DNEQPvty");

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQPvty::Euler3DNEQPvty(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler3DVarSet>(term),
  _library(CFNULL),
  _dhe(3),
  _ye()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  vector<std::string> names(5 + nbSpecies);
  names[0] = "p";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "T";

  // Names for the species
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[5 + ie] = "y" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQPvty::~Euler3DNEQPvty()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQPvty::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQPvty::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DNEQPvty::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQPvty::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::computePhysicalData(const State& state,
					 RealVector& data)
{
  CFreal p = state[0];
  CFreal T = state[4];

  const RealVector& refData = getModel()->getReferencePhysicalData();

  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*getModel()->getTempRef();

  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // Set the species
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = state[5 + ie];
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
    const CFreal w = state[3];
    const CFreal V2 = u*u + v*v + w*w;
    
    data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
    data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2; 
    _library->frozenGammaAndSoundSpeed(Tdim, pdim, rhoDim, 
				       data[EulerTerm::GAMMA], 
				       data[EulerTerm::A], CFNULL);
    data[EulerTerm::V] = sqrt(V2);
    data[EulerTerm::VX] = u;
    data[EulerTerm::VY] = v;
    data[EulerTerm::VZ] = w;
  }
  
  if (iEqSS == 1) {
    // here we reuse the previous composition: we neglect the
    // effect of Ye perturbation on the density rho
    data[EulerTerm::RHO] = _library->density(Tdim,pdim)/
      refData[EulerTerm::RHO];

    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal w = state[3];
    const CFreal V2 = u*u + v*v + w*w;
    data[EulerTerm::V] = sqrt(V2);
    data[EulerTerm::VX] = u;
    data[EulerTerm::VY] = v;
    data[EulerTerm::VZ] = w;
  }

  const CFuint nbEqs = eqSS.getNbEqsSS();
  if ((nbEqs == 5 + nbSpecies) || (iEqSS == 1)) {
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      data[firstSpecies + ie] = state[5 + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::computeStateFromPhysicalData(const RealVector& data,
						 State& state)
{
  state[0] = data[EulerTerm::P];
  state[1] = data[EulerTerm::VX];
  state[2] = data[EulerTerm::VY];
  state[3] = data[EulerTerm::VZ];
  state[4] = data[EulerTerm::T];

  // Set the species
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    state[5 + ie] = data[firstSpecies + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DNEQPvty::getSpeed(const State& state) const
{
  return sqrt(state[1]*state[1] + state[2]*state[2] + state[3]*state[3]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::P];
  result[1] = state[1]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::V];
  result[4] = state[4]*getModel()->getTempRef();

  // Set the species
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    result[5 + ie] = state[5 + ie];
  }

}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::setAdimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::P];
  result[1] = state[1]/refData[EulerTerm::V];
  result[2] = state[2]/refData[EulerTerm::V];
  result[3] = state[3]/refData[EulerTerm::V];
  result[4] = state[4]/(getModel()->getTempRef());

  // Set the species
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    result[5 + ie] = state[5 + ie];
  }

}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  const CFreal u = state[1]*refData[EulerTerm::V];
  const CFreal v = state[2]*refData[EulerTerm::V];
  const CFreal w = state[3]*refData[EulerTerm::V];

  CFreal p = state[0]*refData[EulerTerm::P];
  CFreal T = state[4]*getModel()->getTempRef();

  result[0] = p;
  result[1] = u;
  result[2] = v;
  result[3] = w;
  result[4] = T;

  // Set the species
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = result[5 + ie] = state[5 + ie];
  }

  _library->setSpeciesFractions(_ye);

  // density, total enthalpy, Mach, Xi*nbSpecies
  extra.resize(3);

  const CFreal V2 = u*u + v*v + w*w;
  _library->setDensityEnthalpyEnergy(T, p,_dhe);
  extra[0] = _dhe[0];
  extra[1] = _dhe[1] + 0.5*V2;
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  _library->frozenGammaAndSoundSpeed(T,p,extra[0], gamma, a, CFNULL);
  extra[2] = sqrt(V2)/a;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler3DNEQPvty::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());

  vector<std::string> names(3);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::setup()
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

void Euler3DNEQPvty::computePerturbedPhysicalData
(const State& state, const RealVector& bData, 
 RealVector& data, CFuint iVar) 
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());
  
  if (iVar == 0 || iVar >= 4) { // p + dp // T + dT // Y + dY
    computePhysicalData(state, data);
  }
  else if(iVar == 1 || iVar == 2 || iVar == 3) { // u + du
    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal w = state[3];
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
    
    // Terms for species
    const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
    const CFuint nbSpecies = getModel()->getNbScalarVars(0);
    for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
      data[firstSpecies + ie]= bData[firstSpecies + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvty::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  velIDs.resize(3); 
  velIDs[XX] = 1; velIDs[YY] = 2;  velIDs[ZZ] = 3;      
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
