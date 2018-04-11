#include "NEQ.hh"
#include "Euler2DNEQPivtTvToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQPivtTvToCons, VarSetTransformer, NEQModule, 1>
euler2DNEQPivtTvToConsProvider("Euler2DNEQPivtTvToCons");

Environment::ObjectProvider<Euler2DNEQPivtTvToCons, VarSetTransformer, NEQModule, 1>
euler2DNEQPivtToConsProvider("Euler2DNEQPivtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPivtTvToCons::Euler2DNEQPivtTvToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _dhe(),
  _tvDim(),
  _Rspecies()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPivtTvToCons::~Euler2DNEQPivtTvToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTvToCons::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const RealVector& refData =  _model->getReferencePhysicalData();
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  
  // unused // const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  // Set the mixture density (sum of the partial densities)
  CFreal T = state[nbSpecies + 2];
  CFreal Te = T;
  CFreal Tdim = T*refData[EulerTerm::T];
  
  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 3; 
  if (nbTv > 0) {
    _tvDim.resize(nbTv);
    for (CFuint ie = 0; ie < nbTv; ++ie) {
      _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
    }
    Te = library->getTe(Tdim, &_tvDim[0])/refData[EulerTerm::T];
  }

  _Rspecies.resize(nbSpecies); 
  library->setRiGas(_Rspecies);
  
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal Ti = (ie > 0) ? T : Te;
    result[ie] = state[ie]/(_Rspecies[ie]*Ti);
    rho += result[ie];
  }
  
  // Set the species
  _ye.resize(nbSpecies);
  const CFreal ovRho = 1./rho;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = result[ie]*ovRho;
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
    
  const CFreal u = state[nbSpecies];
  const CFreal v = state[nbSpecies + 1];
  const CFreal V2 = u*u + v*v;
  _dhe.resize(3 + nbTv);
  
  CFreal p = 0.0;
  for (CFuint ie = 0;  ie < nbSpecies; ++ie) {
    p += state[ie]*refData[EulerTerm::P];
  }
  
  const CFreal ovHref = 1./refData[EulerTerm::H];
  
  // vibrational temperatures
  CFreal pdim = _model->getPressureFromState(p);
  if (nbTv == 0) {
    library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
  }
  else {
    library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim,_dhe);
    const CFuint nbTe = library->getNbTe();
    const CFuint nbTvH = nbTv - nbTe; 
    cf_assert(nbTvH >= 0);
    
    // data stores the molecular vibrational energy multiplied 
    // by the molecules mass fractions
    if (nbTvH > 0) {
      CFLog(DEBUG_MAX, "Euler2DNEQPivtTvToCons::transform() => nbTvH > 0\n"); 
      for(CFuint i = 0; i < nbTvH; ++i) {
	result[startTv + i] = rho*_dhe[3 + i]*ovHref; 
      }
    }
    
    if (nbTe == 1) {
      CFLog(DEBUG_MAX, "Euler2DNEQPivtTvToCons::transform() => Te == 1\n"); 
      result[startTv + nbTvH] = rho*_dhe[3 + nbTvH]*ovHref;
    }
  }
  
  result[nbSpecies] = rho*u;
  result[nbSpecies + 1] = rho*v;
  result[nbSpecies + 2] = rho*(_dhe[2]*ovHref + 0.5*V2);
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPivtTvToCons::transformFromRef(const RealVector& data, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFreal rho = data[EulerTerm::RHO];

  //   const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  //   const RealVector& refData =  _model->getReferencePhysicalData();
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*data[firstSpecies + ie];
  }
  
  result[nbSpecies] = rho*data[EulerTerm::VX];
  result[nbSpecies + 1] = rho*data[EulerTerm::VY];
  result[nbSpecies + 2] = rho*data[EulerTerm::E];
  
  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  if (nbTv > 0) {
    const CFuint startTv = nbSpecies + 3;
    const CFuint firstTv = _model->getFirstScalarVar(1);
    for(CFuint ie = 0; ie < nbTv; ++ie) {
      // vibrational temperatures (in data there is y_i*E_i)
      result[startTv + ie] = rho*data[firstTv + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
