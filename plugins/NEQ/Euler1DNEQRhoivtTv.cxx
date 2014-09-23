#include "NEQ.hh"
#include "Euler1DNEQRhoivtTv.hh"
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

Environment::ObjectProvider<Euler1DNEQRhoivtTv, ConvectiveVarSet, NEQModule, 1>
euler1DNEQRhoivtTvProvider("Euler1DNEQRhoivtTv");

//////////////////////////////////////////////////////////////////////////////
      
Euler1DNEQRhoivtTv::Euler1DNEQRhoivtTv(Common::SafePtr<BaseTerm> term) :
  Euler1DNEQRhoivt(term),
  _tvDim(),
  _moleculesIDs()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  _tvDim.resize(nbTv);
 
  vector<std::string> names(nbSpecies + 2 + nbTv);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "rho" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "u";
  names[nbSpecies + 1] = "T";
  
  const CFuint startTv = nbSpecies + 2; 
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    names[startTv + ie] = "Tv" + StringOps::to_str(ie);
  }
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler1DNEQRhoivtTv::~Euler1DNEQRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::computeEigenValuesVectors(RealMatrix& rightTv,
					       RealMatrix& leftTv,
					       RealVector& eValues,
					       const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler1DNEQRhoivtTv::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler1DNEQRhoivtTv::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::splitJacobian(RealMatrix& jacobPlus,
				    RealMatrix& jacobMin,
				    RealVector& eValues,
				    const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler1DNEQRhoivtTv::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::setThermodynamics(CFreal rho, 
					   const State& state, 
					   RealVector& data)
{ 
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const RealVector& refData = getModel()->getReferencePhysicalData();
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal T = state[getTempID(nbSpecies)];
  CFreal Tdim = T*refData[EulerTerm::T];
  
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFuint startTv = nbSpecies + 2;
  
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
  }
  
  CFreal pdim = _library->pressure(rhodim, Tdim, &_tvDim[0]);
  CFreal p = pdim/refData[EulerTerm::P];
  
  // unused //  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  data[EulerTerm::P] = p;
  data[EulerTerm::T] = T;
  data[EulerTerm::RHO] = rho;
  
  if (!_skipEnergyData) {
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

}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::computeStateFromPhysicalData(const RealVector& data,
						 State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler1DNEQRhoivtTv::computeStateFromPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::setDimensionalValues(const State& state,
					      RealVector& result)
{
  // first call the parent
  Euler1DNEQRhoivt::setDimensionalValues(state,result);
  
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 2;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    result[startTv + ie] = state[startTv + ie]*refData[EulerTerm::T];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::setAdimensionalValues(const Framework::State& state,
					       RealVector& result)
{
  // first call the parent
  Euler1DNEQRhoivt::setAdimensionalValues(state,result);
  
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 2;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    result[startTv + ie] = state[startTv + ie]/refData[EulerTerm::T];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  Euler1DNEQRhoivtTv::setDimensionalValues(state,result);
      
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
  CFreal Tdim = state[nbSpecies + 1]*refData[EulerTerm::T];
  // set the vibrational temperature
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 2;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    result[startTv + ie] = state[startTv + ie]*refData[EulerTerm::T];
    _tvDim[ie] = result[startTv + ie];
  }
    
  CFreal pdim = _library->pressure(rhodim, Tdim, &_tvDim[0]);
  _library->setDensityEnthalpyEnergy(Tdim,_tvDim,pdim, _dhe);
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  _library->frozenGammaAndSoundSpeed(Tdim,pdim,rhodim, gamma, a, &_tvDim);
    
  extra.resize(4);
  extra[0] = rho*refData[EulerTerm::RHO];
  extra[1] = _dhe[1] + 0.5*V2;
  extra[2] = sqrt(V2)/a;
  extra[3] = pdim;
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTv::setup()
{
  Euler1DNEQRhoivt::setup();
  
  // set the IDs for the molecules
  _library->setMoleculesIDs(_moleculesIDs);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
