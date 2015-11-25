#include "NEQ.hh"
#include "Euler3DNEQRhoivt.hh"
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

Environment::ObjectProvider<Euler3DNEQRhoivt, ConvectiveVarSet, NEQModule, 1>
euler3DNEQRhoivtProvider("Euler3DNEQRhoivt");

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQRhoivt::Euler3DNEQRhoivt(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler3DVarSet>(term),
  _library(CFNULL),
  _dhe(),
  _ye(),
  _Rspecies()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  vector<std::string> names(nbSpecies + 4);

  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "rho" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "u";
  names[nbSpecies + 1] = "v";
  names[nbSpecies + 2] = "w";
  names[nbSpecies + 3] = "T";

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQRhoivt::~Euler3DNEQRhoivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::computeEigenValuesVectors(RealMatrix& rightEv,
						 RealMatrix& leftEv,
						 RealVector& eValues,
						 const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQRhoivt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DNEQRhoivt::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQRhoivt::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQRhoivt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::computePhysicalData(const State& state, RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie];
  }

  // Set the species
  const CFreal ovRho = 1./rho;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = state[ie]*ovRho;
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
  const CFreal w = state[nbSpecies + 2];
  const CFreal V2 = u*u + v*v +w*w;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::VZ] = w;

  setThermodynamics(rho, state, data);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::computeStateFromPhysicalData(const RealVector& data,
					   State& state)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    state[ie] = data[EulerTerm::RHO]*data[firstSpecies + ie];
  }

  state[nbSpecies]     = data[EulerTerm::VX];
  state[nbSpecies + 1] = data[EulerTerm::VY];
  state[nbSpecies + 2] = data[EulerTerm::VZ];
  state[nbSpecies + 3] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DNEQRhoivt::getSpeed(const State& state) const
{
  const CFuint uID = getModel()->getNbScalarVars(0);
  const CFuint vID = uID + 1;
  const CFuint wID = uID + 2;
  return sqrt(state[uID]*state[uID] + state[vID]*state[vID] + state[wID]*state[wID]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::setDimensionalValues(const State& state,
					    RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]*refData[EulerTerm::RHO];;
  }

  result[nbSpecies]     = state[nbSpecies]*refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]*refData[EulerTerm::V];
  result[nbSpecies + 2] = state[nbSpecies + 2]*refData[EulerTerm::V];
  result[nbSpecies + 3] = state[nbSpecies + 3]*refData[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::setAdimensionalValues(const Framework::State& state,
						 RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]/refData[EulerTerm::RHO];;
  }

  result[nbSpecies]     = state[nbSpecies]/refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]/refData[EulerTerm::V];
  result[nbSpecies + 2] = state[nbSpecies + 2]/refData[EulerTerm::V];
  result[nbSpecies + 3] = state[nbSpecies + 3]/refData[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  // first set the state variables
  Euler3DNEQRhoivt::setDimensionalValues(state,result);

  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie]; // state[ie] is partial density
    _ye[ie] = state[ie];
  }
  _ye /= rho;

  // set the current species fractions in the thermodynamic library
  // this has to be done before computing any other thermodynamic quantity !!!
  _library->setSpeciesFractions(_ye);
  
  CFreal Tdim = state[getTempID(nbSpecies)]*refData[EulerTerm::T];
  CFreal* rhoi = &const_cast<State&>(state)[0];
  _library->setState(rhoi, &Tdim); 
  
  const CFreal u = result[nbSpecies];
  const CFreal v = result[nbSpecies+1];
  const CFreal w = result[nbSpecies+2];
  const CFreal V2 = u*u + v*v + w*w;
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal pdim = _library->pressure(rhodim, Tdim, CFNULL);
  
  _library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  _library->frozenGammaAndSoundSpeed(Tdim,pdim,rhodim, gamma, a, CFNULL);
  
  extra.resize(5);
  extra[0] = rhodim;
  extra[1] = _dhe[1] + 0.5*V2;
  extra[2] = sqrt(V2)/a;
  extra[3] = pdim; 
  extra[4] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler3DNEQRhoivt::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());

  vector<std::string> names(5);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";
  names[3] = "p";
  names[4] = "gamma";
  
  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::setup()
{
  MultiScalarVarSet<Euler3DVarSet>::setup();

  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // set the equation set data for each of the equation subsets
  // first equation subset
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData().resize(2);
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData()[0].setup(0,0,nbSpecies);

  // second equation subset
  Euler3DVarSet::getEqSetData().resize(1);
  Euler3DVarSet::getEqSetData()[0].setup(1,nbSpecies,4);

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());

  const CFuint nbTv = getModel()->getNbScalarVars(1);
  // third equation subset
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData()[1].setup(2,nbSpecies + 4,nbTv);

  _dhe.resize(3 + nbTv);
  _ye.resize(nbSpecies);
  
  _Rspecies.resize(nbSpecies);
  _library->setRiGas(_Rspecies);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::computePerturbedPhysicalData
(const State& state, const RealVector& bData, 
 RealVector& data, CFuint iVar) 
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());
  
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = getTempID(nbSpecies);
  if (iVar < nbSpecies || iVar >= TID) { // rhoi + drhoi // T_i + dT_i
    computePhysicalData(state, data);
  }
  else { // u + du // v + dv // w + dw
    const CFreal u = state[nbSpecies];
    const CFreal v = state[nbSpecies + 1];
    const CFreal w = state[nbSpecies + 2];
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
    
    // Terms for the thermal NEQ
    const CFuint nbTv = getModel()->getNbScalarVars(1);
    const CFuint firstTv = getModel()->getFirstScalarVar(1);
    for (CFuint ie = 0; ie < nbTv; ++ie) {
      data[firstTv + ie] = bData[firstTv + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::setThermodynamics(CFreal rho,
					 const State& state,
					 RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const RealVector& refData = getModel()->getReferencePhysicalData();
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal T = state[getTempID(nbSpecies)];
  CFreal Tdim = T*refData[EulerTerm::T];
  CFreal* rhoi = &const_cast<State&>(state)[0];
  _library->setState(rhoi, &Tdim);
  
  CFreal pdim = _library->pressure(rhodim, Tdim, CFNULL);
  const CFreal p = (pdim - getModel()->getPressInf())/refData[EulerTerm::P];
  
  // unused //  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  data[EulerTerm::P] = p; // dp in the case of incompressible flow
  data[EulerTerm::T] = T;
  data[EulerTerm::RHO] = rho;

  if (!_skipEnergyData) {
    _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
    _library->frozenGammaAndSoundSpeed(Tdim, pdim, rhodim,
				 data[EulerTerm::GAMMA],
				 data[EulerTerm::A], CFNULL);
    
    const CFreal V2 = data[EulerTerm::V]*data[EulerTerm::V];
    data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
    data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::computePressureDerivatives(const Framework::State& state, 
						  RealVector& dp)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 3;
  
  dp[TID] = 0.0;
  for (CFuint  i = 0; i < nbSpecies; ++i) {
    dp[i] = _Rspecies[i]*state[TID];
    dp[TID] += _Rspecies[i]*state[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivt::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  const CFuint nbSpecies = _library->getNbSpecies();
  velIDs.resize(3); 
  velIDs[XX] = nbSpecies; velIDs[YY] = nbSpecies + 1;  velIDs[ZZ] = nbSpecies + 2;   
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
