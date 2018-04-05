#include "NEQ.hh"
#include "Euler2DNEQPivt.hh"
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

Environment::ObjectProvider<Euler2DNEQPivt, ConvectiveVarSet, NEQModule, 1>
euler2DNEQPivtProvider("Euler2DNEQPivt");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPivt::Euler2DNEQPivt(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term),
  _library(CFNULL),
  _dhe(),
  _ye(),
  _rhoi(),
  _Rspecies()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  vector<std::string> names(nbSpecies + 3);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "p" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "u";
  names[nbSpecies + 1] = "v";
  names[nbSpecies + 2] = "T";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPivt::~Euler2DNEQPivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::computeEigenValuesVectors(RealMatrix& rightEv,
					       RealMatrix& leftEv,
					       RealVector& eValues,
					       const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPivt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQPivt::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPivt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////
        
void Euler2DNEQPivt::computePhysicalData(const State& state, RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  // Set the mixture density (sum of the partial densities)
  const CFreal T = state[getTempID(nbSpecies)];
  const CFreal Te = getTe(state);
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Ti = (ie > 0) ? T : Te;
    _rhoi[ie] = state[ie]/(_Rspecies[ie]*Ti);
    rho += _rhoi[ie];
  }
  
  // Set the species
  const CFreal ovRho = 1./rho;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = _rhoi[ie]*ovRho;
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  _library->setSpeciesFractions(_ye);
  
  // set the species mass fractions
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    data[firstSpecies + ie] = _ye[ie];
  }
  
  const CFreal u = state[nbSpecies];
  const CFreal v = state[nbSpecies + 1];
  const CFreal V2 = u*u + v*v;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  
  setThermodynamics(rho, state, data);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::computeStateFromPhysicalData(const RealVector& data,
						  State& state)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    state[ie] = data[EulerTerm::RHO]*data[firstSpecies + ie]*
      _Rspecies[ie]*data[EulerTerm::T];
  } 
  
  state[nbSpecies]     = data[EulerTerm::VX];
  state[nbSpecies + 1] = data[EulerTerm::VY];
  state[nbSpecies + 2] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DNEQPivt::getSpeed(const State& state) const
{
  const CFuint uID = getModel()->getNbScalarVars(0);
  const CFuint vID = uID + 1;
  return sqrt(state[uID]*state[uID] + state[vID]*state[vID]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::setDimensionalValues(const State& state,
					    RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]*refData[EulerTerm::P];
  }
  
  result[nbSpecies]     = state[nbSpecies]*refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]*refData[EulerTerm::V];
  result[nbSpecies + 2] = state[nbSpecies + 2]*refData[EulerTerm::T];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::setAdimensionalValues(const Framework::State& state,
					   RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]/refData[EulerTerm::P]; 
  }
  
  result[nbSpecies]     = state[nbSpecies]/refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]/refData[EulerTerm::V];
  result[nbSpecies + 2] = state[nbSpecies + 2]/refData[EulerTerm::T];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  // first set the state variables
  Euler2DNEQPivt::setDimensionalValues(state,result);
  
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  // Set the mixture density (sum of the partial densities)
  CFreal T = state[getTempID(nbSpecies)];
  const CFreal Te = getTe(state);
  
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Ti = (ie > 0) ? T : Te;
    _rhoi[ie] = state[ie]/(_Rspecies[ie]*Ti);
    rho += _rhoi[ie];
  }

  const CFuint TID = getTempID(nbSpecies);
  CFreal* Tvec = &const_cast<State&>(state)[TID];
  _library->setState(&_rhoi[0], &T);
  
  // Set the species
  const CFreal ovRho = 1./rho;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = _rhoi[ie]*ovRho;
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  _library->setSpeciesFractions(_ye);
  
  const CFreal u = result[nbSpecies];
  const CFreal v = result[nbSpecies+1];
  const CFreal V2 = u*u + v*v;
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal Tdim = state[getTempID(nbSpecies)]*refData[EulerTerm::T];
  CFreal p = 0.0;
  for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
    p += state[ie]*refData[EulerTerm::P];
  }
  
  CFreal pdim = getModel()->getPressureFromState(p);
  _library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  _library->frozenGammaAndSoundSpeed(Tdim,pdim,rhodim, gamma, a, CFNULL);
  
  extra.resize(4);
  extra[0] = rhodim;
  extra[1] = _dhe[1] + 0.5*V2;
  extra[2] = sqrt(V2)/a;
  extra[3] = pdim;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DNEQPivt::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());
  
  vector<std::string> names(4);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";
  names[3] = "p";
  
  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::setup()
{
  MultiScalarVarSet<Euler2DVarSet>::setup();
  
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // set the equation set data for each of the equation subsets
  // first equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData().resize(2);
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].setup(0,0,nbSpecies);
  
  // second equation subset
  Euler2DVarSet::getEqSetData().resize(1);
  Euler2DVarSet::getEqSetData()[0].setup(1,nbSpecies,3);
  
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  // third equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[1].setup(2,nbSpecies + 3,nbTv);
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
  
  _dhe.resize(3 + nbTv);
  _ye.resize(nbSpecies);
  _rhoi.resize(nbSpecies);
  _Rspecies.resize(nbSpecies); 
  _library->setRiGas(_Rspecies);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::computePerturbedPhysicalData(const Framework::State& state,
						  const RealVector& bData,
						  RealVector& data,
						  CFuint iVar)
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());
  
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = getTempID(nbSpecies);
  if (iVar < nbSpecies || iVar >= TID) { // rhoi + drhoi // T_i + dT_i
    computePhysicalData(state, data);
  }
  else { // u + du // v + dv
    const CFreal u = state[nbSpecies];
    const CFreal v = state[nbSpecies + 1];
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
    
    // Terms for the thermal NEQ
    const CFuint nbTv = getModel()->getNbScalarVars(1);
    const CFuint firstTv = getModel()->getFirstScalarVar(1);
    for (CFuint ie = 0; ie < nbTv; ++ie) {
	data[firstTv + ie] = bData[firstTv + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::setThermodynamics(CFreal rho, 
				       const State& state, 
				       RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const RealVector& refData = getModel()->getReferencePhysicalData();
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal T = state[getTempID(nbSpecies)];
  CFreal Tdim = T*refData[EulerTerm::T];
  
  CFreal p = 0.0;
  for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
    p += state[ie];
  }
  
  _library->setState(&_rhoi[0], &Tdim);
  
  // unused //  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  data[EulerTerm::P] = p/refData[EulerTerm::P];
  data[EulerTerm::T] = T;
  data[EulerTerm::RHO] = rho;
  
  if (!_skipEnergyData) {
    CFreal pdim = getModel()->getPressureFromState(p);
    _library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
    _library->frozenGammaAndSoundSpeed(Tdim, pdim, rhodim,
				       data[EulerTerm::GAMMA],
				       data[EulerTerm::A], CFNULL);
    
    const CFreal V2 = data[EulerTerm::V]*data[EulerTerm::V];
    const RealVector& refData = getModel()->getReferencePhysicalData();
    data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
    data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivt::computePressureDerivatives(const Framework::State& state, 
						RealVector& dp)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;
  
  dp[TID] = 0.0;
  for (CFuint  i = 0; i < nbSpecies; ++i) {
    dp[i] = 1; // dp/dp_i
    // dp[TID] += _Rspecies[i]*state[i];
  }
  
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPivt::computePressureDerivatives()");
}
      
//////////////////////////////////////////////////////////////////////////////

bool Euler2DNEQPivt::isValid(const RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint i = 0; i < nbSpecies; ++i) {
    if (data[i] < 0.) {
      CFLog(VERBOSE, "Euler2DNEQPivtTv::isValid() => p_" << i << " = " << data[i] << " < 0 !\n");
      return false;
    }
  }
  
  const CFreal T  = data[nbSpecies+2];
  if (T < 1e-8) {
    CFLog(VERBOSE, "Euler2DNEQPivtTv::isValid() => T = " << T << " < 0!\n");
    return false;
  }
  
  return true;
}

//////////////////////////////////////////////////////////////////////////////
      
CFreal Euler2DNEQPivt::getTe(const Framework::State& state)
{
  const CFreal T = state[getTempID(_library->getNbSpecies())];
  return _library->getTe(T, CFNULL);
}

//////////////////////////////////////////////////////////////////////////////
      
void Euler2DNEQPivt::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  const CFuint nbSpecies = _library->getNbSpecies();
  velIDs.resize(2); velIDs[XX] = nbSpecies; velIDs[YY] = nbSpecies + 1; 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
