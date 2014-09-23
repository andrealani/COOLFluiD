#include "NEQ.hh"
#include "Euler2DNEQRhovtXiToPivtTv.hh"
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

Environment::ObjectProvider<Euler2DNEQRhovtXiToPivtTv, VarSetTransformer, NEQModule, 1>
euler2DNEQRhovtXiToPivtTvProvider("Euler2DNEQRhovtXiToPivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhovtXiToPivtTv::Euler2DNEQRhovtXiToPivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _xe(),
  _tvDim(),
  _Rspecies()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhovtXiToPivtTv::~Euler2DNEQRhovtXiToPivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhovtXiToPivtTv::transform(const State& state, State& result)
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
  
  _Rspecies.resize(nbSpecies); 
  library->setRiGas(_Rspecies);
  
  const RealVector& refData =  _model->getReferencePhysicalData();   
  const CFreal T = state[1];
  CFreal Tdim = T*refData[EulerTerm::T];
  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  _tvDim.resize(nbTv);
  const CFuint startTv = nbSpecies + 3; 
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
    result[startTv + ie] = T; 
  }
  
  const CFreal Te = library->getTe(Tdim, &_tvDim[0])/refData[EulerTerm::T];
  const CFreal rho = state[0];
  
  // partial pressures
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    const CFreal Ti = (ie > 0) ? T : Te;
    result[ie] = rho*_ye[ie]*_Rspecies[ie]*Ti;
  }
  
  result[nbSpecies]     = state[2]; // u 
  result[nbSpecies + 1] = state[3]; // v
  result[nbSpecies + 2] = T;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
