#include "NEQKOmega.hh"
#include "Euler2DNEQKOmegaRhoivtToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQKOmegaRhoivtToCons, VarSetTransformer, NEQKOmegaModule, 1>
euler2DNEQKOmegaRhoivtToConsProvider("Euler2DNEQKOmegaRhoivtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQKOmegaRhoivtToCons::Euler2DNEQKOmegaRhoivtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerNEQKOmegaTerm>()),
  _dhe(3),
  _ye()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQKOmegaRhoivtToCons::~Euler2DNEQKOmegaRhoivtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQKOmegaRhoivtToCons::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbSpecies = _model->getNbScalarVars(0);
  //   const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  // unused //  const CFuint iEqSS = eqSS.getEqSS();
  // unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  // partial densities
  _ye.resize(nbSpecies);
  CFreal rho = 0.0;
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = state[ie];
    rho += state[ie];
    _ye[ie] = state[ie];
  }
  _ye /= rho;

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!!
  library->setSpeciesFractions(_ye);

  const RealVector& refData =  _model->getReferencePhysicalData();
  CFreal rhoDim = rho*refData[EulerTerm::RHO];
  CFreal T = state[nbSpecies + 2];
  CFreal Tdim = T*refData[EulerTerm::T];
  CFreal pdim = library->pressure(rhoDim, Tdim, CFNULL);
  const CFreal u = state[nbSpecies];
  const CFreal v = state[nbSpecies + 1];
  const CFreal V2 = u*u + v*v;
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  result[nbSpecies] = rho*u;
  result[nbSpecies + 1] = rho*v;
  result[nbSpecies + 2] = rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
  
  result[nbSpecies + 3] = rho*state[nbSpecies + 3]; // k
  result[nbSpecies + 4] = rho*state[nbSpecies + 4]; // omega
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQKOmegaRhoivtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFreal rho = data[EulerTerm::RHO];
    
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*data[firstSpecies + ie];
  }
  
  result[nbSpecies] = rho*data[EulerTerm::VX];
  result[nbSpecies + 1] = rho*data[EulerTerm::VY];
  result[nbSpecies + 2] = rho*data[EulerTerm::E];

  const CFuint kID = _model->getFirstScalarVar(1);
  result[nbSpecies + 3] = rho*data[kID];
  result[nbSpecies + 4] = rho*data[kID + 1];
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQKOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
