#include "NEQ.hh"
#include "Euler3DNEQRhovtXiToRhoivtTv.hh"
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

Environment::ObjectProvider<Euler3DNEQRhovtXiToRhoivtTv, VarSetTransformer, NEQModule, 1>
Euler3DNEQRhovtXiToRhoivtTvProvider("Euler3DNEQRhovtXiToRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQRhovtXiToRhoivtTv::Euler3DNEQRhovtXiToRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQRhovtXiToRhoivtTv::~Euler3DNEQRhovtXiToRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhovtXiToRhoivtTv::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  _ye.resize(nbSpecies);
  _xe.resize(nbSpecies);
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _xe[ie] = state[5 + ie];
  }
    
  library->getSpeciesMassFractions(_xe, _ye);
    
  const CFreal rho = state[0];
  const CFreal T = state[1];
    
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*_ye[ie];
   // CFLog(INFO, "result[ie] = " <<result[ie] <<"\n");
  }
    
  result[nbSpecies]     = state[2]; // u 
  result[nbSpecies + 1] = state[3]; // v
  result[nbSpecies + 2] = state[4]; // w
  result[nbSpecies + 3] = T;
  
  // set the vibrational temperature
  const CFuint startTv = nbSpecies + 4;
  const CFuint nbTv = _model->getNbScalarVars(1);
  for (CFuint i = 0; i < nbTv; ++i){
    result[startTv + i] = T;
  }   
 // CFLog(INFO, "Euler3DNEQRhovtXiToRhoivtTv" <<"\n");
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
