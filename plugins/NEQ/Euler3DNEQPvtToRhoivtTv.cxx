#include "NEQ.hh"
#include "Euler3DNEQPvtToRhoivtTv.hh"
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

Environment::ObjectProvider<Euler3DNEQPvtToRhoivtTv, VarSetTransformer, NEQModule, 1>
euler3DNEQPvtToRhoivtTvProvider("Euler3DNEQPvtToRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQPvtToRhoivtTv::Euler3DNEQPvtToRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe(),
  _masses()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQPvtToRhoivtTv::~Euler3DNEQPvtToRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQPvtToRhoivtTv::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
 
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  _ye.resize(nbSpecies);
  _xe.resize(nbSpecies);
  _masses.resize(nbSpecies);
  library->getMolarMasses(_masses);

  CFdouble T = state[4];
  CFdouble p = state[0];
  CFreal mm = 0;

  library->setComposition(T, p, &_xe); // to clarify this part for the computation of _xe

  for (CFuint i = 0; i < nbSpecies; i++)
    {
      mm += _xe[i]*_masses[i];
    }
 
  const CFreal Rmix = library->getRgas()/mm;
  const CFreal rho = state[0]/(Rmix*T); 

  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = _xe[ie]*_masses[ie]/mm; 
    result[ie] = rho*_ye[ie];
  }
   
  // TRY THIS !!!!! 
  /*_ye[0] = _ye[1] = _ye[3] = 0.0; 
  _ye[2] = 0.767; _ye[4] = 0.233;
   
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*_ye[ie];
   } */

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  cf_assert(result.size() == state.size());
   
  result[nbSpecies]     = state[1];
  result[nbSpecies + 1] = state[2];
  result[nbSpecies + 2] = state[3];
  result[nbSpecies + 3] = T;
  
  // set the vibrational temperature
  const CFuint startTv = nbSpecies + 4;
  const CFuint nbTv = _model->getNbScalarVars(1);
  for (CFuint i = 0; i < nbTv; ++i){
     result[startTv + i] = (T>2000.) ? 2000. : T;
     // result[startTv + i] = 248.528;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
