#include "NEQ.hh"
#include "Euler3DNEQPvtToRhoivt.hh"
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

Environment::ObjectProvider<Euler3DNEQPvtToRhoivt, VarSetTransformer, NEQModule, 1>
euler3DNEQPvtToRhoivtProvider("Euler3DNEQPvtToRhoivt");

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQPvtToRhoivt::Euler3DNEQPvtToRhoivt
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQPvtToRhoivt::~Euler3DNEQPvtToRhoivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvtToRhoivt::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
 
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  _ye.resize(nbSpecies);
  if (nbSpecies == 2) {
    // we assume N2 mixture
    _ye[0] = 0.;
    _ye[1] = 1.;
  }
  if (nbSpecies == 5) {
    // we assume air5
    _ye[0] = _ye[1] = _ye[3] = 0.;
    _ye[2] = 0.767082;
    _ye[4] = 0.232917;    
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  
  cf_assert(result.size() == state.size());
  
  const CFreal Rmix = library->getRgas()/library->getMMass();
  const CFreal rho = state[0]/(Rmix*state[4]);
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*_ye[ie];
  }
  
  result[nbSpecies]     = state[1];
  result[nbSpecies + 1] = state[2];
  result[nbSpecies + 2] = state[3];
  result[nbSpecies + 3] = state[4];
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
