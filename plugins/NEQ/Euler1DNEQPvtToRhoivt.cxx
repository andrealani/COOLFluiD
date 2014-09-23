#include "NEQ.hh"
#include "Euler1DNEQPvtToRhoivt.hh"
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

Environment::ObjectProvider<Euler1DNEQPvtToRhoivt, VarSetTransformer, NEQModule, 1>
euler1DNEQPvtToRhoivtProvider("Euler1DNEQPvtToRhoivt");

//////////////////////////////////////////////////////////////////////////////

Euler1DNEQPvtToRhoivt::Euler1DNEQPvtToRhoivt
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler1DNEQPvtToRhoivt::~Euler1DNEQPvtToRhoivt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQPvtToRhoivt::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  _ye.resize(nbSpecies);
  _xe.resize(nbSpecies);
  
  RealVector masses;
  masses.resize(nbSpecies);
  library-> getMolarMasses(masses);
  
  CFdouble T = state[2];
  CFdouble p = state[0];
  CFreal mm = 0;
  
  library->setComposition(T, p, &_xe); 
  
  for (CFuint i = 0; i < nbSpecies; i++)
    {
      mm += _xe[i]*masses[i];
    }
  
  const CFreal Rmix = library->getRgas()/mm;
  const CFreal rho = state[0]/(Rmix*state[2]); 
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = _xe[ie]*masses[ie]/mm; 
    result[ie] = rho*_ye[ie];
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  cf_assert(result.size() == state.size());
  
  result[nbSpecies]     = state[1];
  result[nbSpecies + 1] = state[2];
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
