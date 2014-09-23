#include "ICP/ICPNEQ.hh"
#include "ICPNEQ2DPvtToPivtTv.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ICPNEQ2DPvtToPivtTv, VarSetTransformer, ICPNEQModule,1>
icpNEQ2DPvtToPivtTvProvider("ICPNEQ2DPvtToPivtTv");

//////////////////////////////////////////////////////////////////////

ICPNEQ2DPvtToPivtTv::ICPNEQ2DPvtToPivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe(),
  _masses()
{
}
      
//////////////////////////////////////////////////////////////////////
      
ICPNEQ2DPvtToPivtTv::~ICPNEQ2DPvtToPivtTv()
{
}

//////////////////////////////////////////////////////////////////////

void ICPNEQ2DPvtToPivtTv::transform(const State& state, State& result)
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
  
  CFreal T = state[3];
  const CFreal p = state[0];
  CFreal pdim = _model->getPressureFromState(p);
  CFreal mm = 0.;
  library->setComposition(T, pdim, &_xe);
  for (CFuint i = 0; i < nbSpecies; i++) {
    mm += _xe[i]*_masses[i];
  }
  
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = _xe[ie]*_masses[ie]/mm; 
    // p_i=x_i*dp, dp_i=x_i*dp
    result[ie] = _xe[ie]*p;
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
    
  const CFuint startE = dim+2; 
  result[result.size()-2] = state[startE];
  result[result.size()-1] = state[startE+1];
}

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
