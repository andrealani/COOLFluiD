#include "NEQ.hh"
#include "Euler2DNEQRhovtXiToRhoivtTv.hh"
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

Environment::ObjectProvider<Euler2DNEQRhovtXiToRhoivtTv, VarSetTransformer, NEQModule, 1>
euler2DNEQRhovtXiToRhoivtTvProvider("Euler2DNEQRhovtXiToRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhovtXiToRhoivtTv::Euler2DNEQRhovtXiToRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhovtXiToRhoivtTv::~Euler2DNEQRhovtXiToRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhovtXiToRhoivtTv::transform(const State& state, State& result)
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
    _xe[ie] = state[4 + ie];
  }
    
  library->getSpeciesMassFractions(_xe, _ye);
    
  const CFreal rho = state[0];
  const CFreal T = state[1];
    
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*_ye[ie];
  }
    
  result[nbSpecies]     = state[2]; // u 
  result[nbSpecies + 1] = state[3]; // v
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
