#include "NEQ.hh"
#include "Euler2DNEQPivtTv.hh"
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

Environment::ObjectProvider<Euler2DNEQPivtTv, ConvectiveVarSet, NEQModule, 1>
euler2DNEQPivtTvProvider("Euler2DNEQPivtTv");

//////////////////////////////////////////////////////////////////////////////
      
Euler2DNEQPivtTv::Euler2DNEQPivtTv(Common::SafePtr<BaseTerm> term) :
  Euler2DNEQPivt(term),
  _tvDim(),
  _moleculesIDs()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  _tvDim.resize(nbTv);
  
  vector<std::string> names(nbSpecies + 3 + nbTv);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "p" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "u";
  names[nbSpecies + 1] = "v";
  names[nbSpecies + 2] = "T";
  
  const CFuint startTv = nbSpecies + 3; 
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    names[startTv + ie] = "Tv" + StringOps::to_str(ie);
  }
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPivtTv::~Euler2DNEQPivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::computeEigenValuesVectors(RealMatrix& rightTv,
						 RealMatrix& leftTv,
						 RealVector& eValues,
						 const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPivtTv::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQPivtTv::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::splitJacobian(RealMatrix& jacobPlus,
				     RealMatrix& jacobMin,
				     RealVector& eValues,
				     const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPivtTv::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::setThermodynamics(CFreal rho, 
					 const State& state, 
					 RealVector& data)
{ 
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const RealVector& refData = getModel()->getReferencePhysicalData();
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal T = state[getTempID(nbSpecies)];
  CFreal Tdim = T*refData[EulerTerm::T];

  CFLog(DEBUG_MIN, "Euler2DNEQPivtTv::setThermodynamics() START\n");
  
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFuint startTv = nbSpecies + 3;
 
  CFLog(DEBUG_MIN, "nbTv=" << nbTv << ", ");
  CFLog(DEBUG_MIN, "firstTv=" << firstTv << ", ");
  CFLog(DEBUG_MIN, "startTv=" << startTv << "\n");

  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
     CFLog(DEBUG_MIN, "tvDim[" << ie << "] = " << _tvDim[ie] << " ");
  }
  CFLog(DEBUG_MIN, "tvDim= " << _tvDim << "\n");
 
  CFreal p = 0.0;
  for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
    p += state[ie];
  }

  CFLog(DEBUG_MIN, "rhodim=" << rhodim << ", ");
  CFLog(DEBUG_MIN, "Tdim=" << Tdim << ", ");
  CFLog(DEBUG_MIN, "p=" << p << "\n");
   
  // unused //  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  data[EulerTerm::P] = p/refData[EulerTerm::P];
  data[EulerTerm::T] = T;
  data[EulerTerm::RHO] = rho;
  
  if (!_skipEnergyData) {
    CFreal pdim = getModel()->getPressureFromState(p);
    _library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim,_dhe,_extraData);
    
    const CFuint nbTe = _library->getNbTe();
    const CFuint nbTvH = nbTv - nbTe;

    // data stores the moleculare vibrational energy multiplied 
    // by the molecules mass fractions
    if (nbTvH != 0) {
       for (CFuint ie = 0; ie < nbTvH; ++ie) {
           data[firstTv + ie] = _dhe[3 + ie]/refData[EulerTerm::H]; 
       } 
    }
 
    if (nbTe == 1) {
      data[firstTv + nbTvH] = _dhe[3 + nbTvH]/refData[EulerTerm::H];
    }
    
    _library->frozenGammaAndSoundSpeed(Tdim, pdim, rhodim,
				       data[EulerTerm::GAMMA],
				       data[EulerTerm::A], &_tvDim);
    
    const CFreal V2 = data[EulerTerm::V]*data[EulerTerm::V];
    data[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
    data[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
  }
 
  CFLog(DEBUG_MIN, "data=" << data << "\n");
  CFLog(DEBUG_MIN, "Euler2DNEQPivtTv::setThermodynamics() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::computeStateFromPhysicalData(const RealVector& data,
						 State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQPivtTv::computeStateFromPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::setDimensionalValues(const State& state,
					      RealVector& result)
{
  // first call the parent
  Euler2DNEQPivt::setDimensionalValues(state,result);
  
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 3;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    result[startTv + ie] = state[startTv + ie]*refData[EulerTerm::T];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::setAdimensionalValues(const Framework::State& state,
					       RealVector& result)
{
  // first call the parent
  Euler2DNEQPivt::setAdimensionalValues(state,result);
  
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 3;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    result[startTv + ie] = state[startTv + ie]/refData[EulerTerm::T];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  Euler2DNEQPivtTv::setDimensionalValues(state,result);
      
  const RealVector& refData = getModel()->getReferencePhysicalData();
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
      
  const CFreal u = result[nbSpecies];
  const CFreal v = result[nbSpecies+1];
  const CFreal V2 = u*u + v*v;
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal Tdim = state[nbSpecies + 2]*refData[EulerTerm::T];
  // set the vibrational temperature
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 3;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    result[startTv + ie] = state[startTv + ie]*refData[EulerTerm::T];
    _tvDim[ie] = result[startTv + ie];
  }
  
  CFreal p = 0.0;
  for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
    p += state[ie]*refData[EulerTerm::P];
  }
  
  CFreal pdim = getModel()->getPressureFromState(p);
  _library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim, _dhe);
  CFLog(DEBUG_MIN, "Euler2DNEQPivtTv::setDimensionalValuesPlusExtraValues() START\n");
  CFLog(DEBUG_MIN, "Tdim=" << Tdim << ", _tvDim = " << _tvDim << ", pdim = " << pdim << ", dhe = " << _dhe << "\n");
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  _library->frozenGammaAndSoundSpeed(Tdim,pdim,rhodim, gamma, a, &_tvDim);
  CFLog(DEBUG_MIN, "gamma=" << gamma << ", a= " << a << "\n");
  CFLog(DEBUG_MIN, "Euler2DNEQPivtTv::setDimensionalValuesPlusExtraValues() END\n");

  extra.resize(4);
  extra[0] = rho*refData[EulerTerm::RHO];
  extra[1] = _dhe[1] + 0.5*V2;
  extra[2] = sqrt(V2)/a;
  extra[3] = pdim;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::setup()
{
  Euler2DNEQPivt::setup();
  
  // set the IDs for the molecules
  _library->setMoleculesIDs(_moleculesIDs);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTv::computePressureDerivatives(const Framework::State& state, 
						  RealVector& dp)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;
  
  dp[TID] = 0.0;
  for (CFuint  i = 0; i < nbSpecies; ++i) {
    dp[i] = 1.;
  }
  
  if (_library->presenceElectron()) {
    const CFuint TeID = (_library->getNbTe() == 0) ? 
      TID + 1 + _library->getElectrEnergyID() : state.size() - 1;
    dp[TeID] = 0.;
  }
}

//////////////////////////////////////////////////////////////////////////////

bool Euler2DNEQPivtTv::isValid(const RealVector& data)
{
  bool correct = true;

  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  //These are indices!
  //  const CFuint Ux = nbSpecies;
  // const CFuint Uy = nbSpecies+1;
  const CFuint temp = nbSpecies+2;
  const CFuint tempVib = nbSpecies+3;

  //These are the actual values:
  const CFreal T = data[temp];
  const CFreal Tv = data[tempVib];
  //
  
  RealVector molecularMasses(nbSpecies);
  _library->getMolarMasses(molecularMasses);

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

  CFreal denom = 0.;
  CFreal riovermi  = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    riovermi += data[i]/molecularMasses[i];
    const CFreal yOvM = data[i]/molecularMasses[i];
    denom += yOvM*((Rgas*fCoeff[i]));
  }

  //Compute sound speed:
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal sigmai = data[i]/molecularMasses[i];
    numBeta += sigmai;
    denBeta += sigmai*fCoeff[i];
  }

  const CFreal beta = numBeta/denBeta;
  
  CFreal  rho = 0.;
  // Set the species:
  for (CFuint ie = 0; ie < nbSpecies; ++ie){
    rho += data[ie];
    if( data[ie] < 0. ){
      return correct = false;
    }
  }

  const CFreal p = Rgas*riovermi*T;
  const CFreal a = std::sqrt((1+beta)*p/rho);

  if( (p < 0.) || (T < 0.) || (Tv < 0.) || (a < 0.) ){
  return correct = false;
  }

  return correct;
}

//////////////////////////////////////////////////////////////////////////////
   
CFreal Euler2DNEQPivtTv::getTe(const Framework::State& state)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFreal T = state[nbSpecies + 2];
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 3;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie];
  }
  
  const CFreal Te = _library->getTe(T, &_tvDim[0]);
  CFLog(DEBUG_MIN, "Euler2DNEQPivtTv::getTe() => Te=" << Te << "\n"); 
  cf_assert(Te > 0.);

  return Te;
}

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
