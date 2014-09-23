#include "NEQ.hh"
#include "Euler1DNEQRhoivt.hh"
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

Environment::ObjectProvider<Euler1DNEQRhoivt, ConvectiveVarSet, NEQModule, 1>
euler1DNEQRhoivtProvider("Euler1DNEQRhoivt");

//////////////////////////////////////////////////////////////////////////////

Euler1DNEQRhoivt::Euler1DNEQRhoivt(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler1DVarSet>(term),
  _library(CFNULL),
  _dhe(),
  _ye()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  vector<std::string> names(nbSpecies + 2);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "rho" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "u";
  names[nbSpecies + 1] = "T";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler1DNEQRhoivt::~Euler1DNEQRhoivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::computeEigenValuesVectors(RealMatrix& rightEv,
						 RealMatrix& leftEv,
						 RealVector& eValues,
						 const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler1DNEQRhoivt::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler1DNEQRhoivt::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler1DNEQRhoivt::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::computePhysicalData(const State& state, RealVector& data)
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
  const CFreal V2 = u*u; 
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  
  setThermodynamics(rho, state, data);
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::computeStateFromPhysicalData(const RealVector& data,
						    State& state)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    state[ie] = data[EulerTerm::RHO]*data[firstSpecies + ie];
  } 
  
  state[nbSpecies]     = data[EulerTerm::VX];
  state[nbSpecies + 1] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler1DNEQRhoivt::getSpeed(const State& state) const
{
  const CFuint uID = getModel()->getNbScalarVars(0);
  return sqrt(state[uID]*state[uID]);
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::setDimensionalValues(const State& state,
					    RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]*refData[EulerTerm::RHO];
  }
  
  result[nbSpecies]     = state[nbSpecies]*refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]*refData[EulerTerm::T];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::setAdimensionalValues(const Framework::State& state,
						 RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]/refData[EulerTerm::RHO]; 
  }
  
  result[nbSpecies]     = state[nbSpecies]/refData[EulerTerm::V];
  result[nbSpecies + 1] = state[nbSpecies + 1]/refData[EulerTerm::T];
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  // first set the state variables
  Euler1DNEQRhoivt::setDimensionalValues(state,result);
  
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

  const CFreal u = result[nbSpecies];
  const CFreal V2 = u*u; 
  CFreal rhodim = rho*refData[EulerTerm::RHO]; 
  CFreal Tdim = state[getTempID(nbSpecies)]*refData[EulerTerm::T];
  CFreal pdim = _library->pressure(rhodim, Tdim, CFNULL);

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

vector<std::string> Euler1DNEQRhoivt::getExtraVarNames() const
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

void Euler1DNEQRhoivt::setup()
{
  MultiScalarVarSet<Euler1DVarSet>::setup();
  
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // set the equation set data for each of the equation subsets
  // first equation subset
  MultiScalarVarSet<Euler1DVarSet>::getEqSetData().resize(2);
  MultiScalarVarSet<Euler1DVarSet>::getEqSetData()[0].setup(0,0,nbSpecies);  
  
  // second equation subset
  Euler1DVarSet::getEqSetData().resize(1);
  Euler1DVarSet::getEqSetData()[0].setup(1,nbSpecies,2);
  
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  // third equation subset
  MultiScalarVarSet<Euler1DVarSet>::getEqSetData()[1].setup(2,nbSpecies + 2,nbTv);  
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
  
  _dhe.resize(3 + nbTv);
  _ye.resize(nbSpecies);
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::computePerturbedPhysicalData
(const State& state, const RealVector& pdataBkp, 
 RealVector& pdata, CFuint iVar) 
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());    
  
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = getTempID(nbSpecies);
  if (iVar < nbSpecies || iVar >= TID) { // rhoi + drhoi // T_i + dT_i
    computePhysicalData(state, pdata);
  }
  else { // u + du 
    const CFreal u = state[nbSpecies];
    const CFreal V2 = u*u; 
    
    pdata[EulerTerm::VX] = u;
    pdata[EulerTerm::RHO] = pdataBkp[EulerTerm::RHO];
    pdata[EulerTerm::V] = sqrt(V2);
    pdata[EulerTerm::P] = pdataBkp[EulerTerm::P];
    pdata[EulerTerm::H] = pdataBkp[EulerTerm::H] +
      0.5* (V2 - pdataBkp[EulerTerm::V]*pdataBkp[EulerTerm::V]);
    pdata[EulerTerm::A] = pdataBkp[EulerTerm::A];
    pdata[EulerTerm::GAMMA] = pdataBkp[EulerTerm::GAMMA];
    pdata[EulerTerm::E] = pdataBkp[EulerTerm::E] +
      0.5* (V2 - pdataBkp[EulerTerm::V]*pdataBkp[EulerTerm::V]);
    
    // Terms for species
    const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
    const CFuint nbSpecies = getModel()->getNbScalarVars(0);
    for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
      pdata[firstSpecies + ie]= pdataBkp[firstSpecies + ie];
    }
    
    // Terms for the thermal NEQ
    const CFuint nbTv = getModel()->getNbScalarVars(1);
    const CFuint firstTv = getModel()->getFirstScalarVar(1);
    for (CFuint ie = 0; ie < nbTv; ++ie) {
      pdata[firstTv + ie] = pdataBkp[firstTv + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivt::setThermodynamics(CFreal rho, 
					 const State& state, 
					 RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const RealVector& refData = getModel()->getReferencePhysicalData();
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal T = state[getTempID(nbSpecies)];
  CFreal Tdim = T*refData[EulerTerm::T];
  CFreal pdim = _library->pressure(rhodim, Tdim, CFNULL);
  CFreal p = pdim/refData[EulerTerm::P];
  
  // unused //  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  data[EulerTerm::P] = p;
  data[EulerTerm::T] = T;
  data[EulerTerm::RHO] = rho;
  
  if (!_skipEnergyData) {
    _library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
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

void Euler1DNEQRhoivt::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  const CFuint nbSpecies = _library->getNbSpecies();
  velIDs.resize(1); velIDs[XX] = nbSpecies;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
