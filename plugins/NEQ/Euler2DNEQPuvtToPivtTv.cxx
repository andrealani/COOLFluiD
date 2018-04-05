#include "NEQ.hh"
#include "Euler2DNEQPuvtToPivtTv.hh"
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

Environment::ObjectProvider<Euler2DNEQPuvtToPivtTv, VarSetTransformer, NEQModule, 1>
euler2DNEQPuvtToPivtTvProvider("Euler2DNEQPuvtToPivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPuvtToPivtTv::Euler2DNEQPuvtToPivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe(),
  _masses()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPuvtToPivtTv::~Euler2DNEQPuvtToPivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPuvtToPivtTv::transform(const State& state, State& result)
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

  // the following has to be adapted for enabling adimensional cases  
  CFdouble T = state[3];
  // state[0] can be "p" or "dp" 
  CFdouble pdim = state[0] + std::max(_model->getPressInfComp(), _model->getPressInf());
  cf_assert(pdim > -1.e-9); 

  CFLog(DEBUG_MIN, "Euler2DNEQPuvtToPivtTv::transform() => pdim = " << pdim << "\n");

  library->setComposition(T, pdim, &_xe);
  
  CFreal mm = 0.;
  for (CFuint i = 0; i < nbSpecies; i++) {
    mm += _xe[i]*_masses[i];
  }
  
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = _xe[ie]*_masses[ie]/mm; 
    // p_i=x_i*dp, dp_i=x_i*dp
    result[ie] = _xe[ie]*pdim;
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  cf_assert(result.size() == state.size());
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint i=0; i < dim; ++i) {
    result[nbSpecies+i] = state[1+i];
  }
  result[nbSpecies + dim] = T;
  
  // set the vibrational temperature
  const CFuint startTv = nbSpecies + dim + 1;
  const CFuint nbTv = _model->getNbScalarVars(1);
  for (CFuint i = 0; i < nbTv; ++i) {
    result[startTv + i] = T;
  }
  
  CFLog(DEBUG_MAX, "Euler2DNEQPuvtToPivtTv::transform() => result = " << result << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
