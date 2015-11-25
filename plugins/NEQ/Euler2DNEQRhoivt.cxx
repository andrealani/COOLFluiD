#include "NEQ.hh"
#include "Euler2DNEQRhoivt.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

#include <fstream>

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

Environment::ObjectProvider<Euler2DNEQRhoivt, ConvectiveVarSet, NEQModule, 1>
euler2DNEQRhoivtProvider("Euler2DNEQRhoivt");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivt::Euler2DNEQRhoivt(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term),
  _library(CFNULL),
  _dhe(),
  _ye()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  vector<std::string> names(nbSpecies + 3);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "rho" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "u";
  names[nbSpecies + 1] = "v";
  names[nbSpecies + 2] = "T";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivt::~Euler2DNEQRhoivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::computeEigenValuesVectors(RealMatrix& rightEv,
						 RealMatrix& leftEv,
						 RealVector& eValues,
						 const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQRhoivt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQRhoivt::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQRhoivt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::computePhysicalData(const State& state, RealVector& data)
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
  const CFreal V2 = u*u + v*v;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  
  setThermodynamics(rho, state, data);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::computeStateFromPhysicalData(const RealVector& data,
					   State& state)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    state[ie] = data[EulerTerm::RHO]*data[firstSpecies + ie];
  } 
  
  state[nbSpecies]     = data[EulerTerm::VX];
  state[nbSpecies + 1] = data[EulerTerm::VY];
  state[nbSpecies + 2] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DNEQRhoivt::getSpeed(const State& state) const
{
  const CFuint uID = getModel()->getNbScalarVars(0);
  const CFuint vID = uID + 1;
  return sqrt(state[uID]*state[uID] + state[vID]*state[vID]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::setDimensionalValues(const State& state,
					    RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]*refData[EulerTerm::RHO];
  }
  
  result[nbSpecies]     = state[nbSpecies]*refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]*refData[EulerTerm::V];
  result[nbSpecies + 2] = state[nbSpecies + 2]*refData[EulerTerm::T];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::setAdimensionalValues(const Framework::State& state,
						 RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]/refData[EulerTerm::RHO]; 
  }
  
  result[nbSpecies]     = state[nbSpecies]/refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]/refData[EulerTerm::V];
  result[nbSpecies + 2] = state[nbSpecies + 2]/refData[EulerTerm::T];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  // first set the state variables
  Euler2DNEQRhoivt::setDimensionalValues(state,result);
  
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
  const CFreal V2 = u*u + v*v;
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

vector<std::string> Euler2DNEQRhoivt::getExtraVarNames() const
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

void Euler2DNEQRhoivt::setup()
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
  
  _Rspecies.resize(nbSpecies);
  
  _library->setRiGas(_Rspecies); 
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::computePerturbedPhysicalData(const Framework::State& state,
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

void Euler2DNEQRhoivt::setThermodynamics(CFreal rho, 
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

void Euler2DNEQRhoivt::computePressureDerivatives(const Framework::State& state, 
						  RealVector& dp)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;
  
  dp[TID] = 0.0;
  for (CFuint  i = 0; i < nbSpecies; ++i) {
    dp[i] = _Rspecies[i]*state[TID];
    dp[TID] += _Rspecies[i]*state[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

bool Euler2DNEQRhoivt::isValid(const RealVector& data)
{
  bool correct = true;

  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  //These are indices!
  const CFuint temp = nbSpecies+2;
  //
  const CFreal T = data[temp];
  //

  // needed for beta coefficient
  RealVector molecularMasses(nbSpecies);
  _library->getMolarMasses(molecularMasses);
  cout << endl << molecularMasses << endl;

  vector<CFuint> moleculeIDs;
  _library->setMoleculesIDs(moleculeIDs);
  vector<bool> flag(nbSpecies, false);
  for (CFuint i = 0; i < moleculeIDs.size(); ++i) {
    flag[moleculeIDs[i]] = true;
  }

  RealVector fCoeff(nbSpecies);
  for (CFuint i = 0; i < nbSpecies; ++i) {
    fCoeff[i] = (flag[i]) ? 2.5 : 1.5;
  }
 
  const CFreal Rgas = _library->getRgas();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
  //

  // Compute rho:
  //Think on accessing species like this:
  //   data[firstSpecies + ie] < 0.
  //in that way could be more general.
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += data[ie];
  }

  const CFreal ovRho = 1./rho;
  //

  CFreal denom = 0.;
  CFreal form  = 0.;
  CFreal riovermi  = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    riovermi += data[i]/molecularMasses[i];
    const CFreal yOvM = data[i]/molecularMasses[i];
    denom += yOvM*((Rgas*fCoeff[i]));
    form += data[i]*eData->enthalpyForm[i];
  }

  const CFreal p = T*Rgas*riovermi;

  //Compute sound speed:
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal sigmai = data[i]/molecularMasses[i];
    numBeta += sigmai;
    denBeta += sigmai*fCoeff[i];
  }

  const CFreal beta = numBeta/denBeta;

  const CFreal a = std::sqrt((1+beta)*p*ovRho);

  if( ( p < 0.) || (T < 0.) || (a < 0.) ){
  return correct = false;
  }

  // Check positivity of the species:

  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    if( data[ie] < 0. ){
      return correct = false;
    }
  }

return correct;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivt::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  const CFuint nbSpecies = _library->getNbSpecies();
  velIDs.resize(2); velIDs[XX] = nbSpecies; velIDs[YY] = nbSpecies + 1; 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
