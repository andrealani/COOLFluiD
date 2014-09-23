#include "NEQ.hh"
#include "Euler2DNEQPuvtToRhoivtTv.hh"
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

Environment::ObjectProvider<Euler2DNEQPuvtToRhoivtTv, VarSetTransformer, NEQModule, 1>
euler2DNEQPuvtToRhoivtTvProvider("Euler2DNEQPuvtToRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPuvtToRhoivtTv::Euler2DNEQPuvtToRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe(),
  _masses()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPuvtToRhoivtTv::~Euler2DNEQPuvtToRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPuvtToRhoivtTv::transform(const State& state, State& result)
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
  
  CFdouble T = state[3];
  // state[0] can be "p" or "dp" 
  CFdouble pdim = state[0] + std::max(_model->getPressInfComp(), _model->getPressInf());
  
  library->setComposition(T, pdim, &_xe);
  
  CFreal mm = 0.;
  for (CFuint i = 0; i < nbSpecies; i++) {
    mm += _xe[i]*_masses[i];
  }
  
  const CFreal Rmix = library->getRgas()/mm;
  const CFreal rho = pdim/(Rmix*T); 
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = _xe[ie]*_masses[ie]/mm; 
    result[ie] = rho*_ye[ie];
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  cf_assert(result.size() == state.size());
   
  result[nbSpecies]     = state[1];
  result[nbSpecies + 1] = state[2];
  result[nbSpecies + 2] = T;
  
  // set the vibrational temperature
  const CFuint startTv = nbSpecies + 3;
  const CFuint nbTv = _model->getNbScalarVars(1);
  for (CFuint i = 0; i < nbTv; ++i){
    result[startTv + i] = T;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
