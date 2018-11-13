#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolumePoisson/PureDiffFlux.hh"
#include "FiniteVolumePoisson/FiniteVolumePoisson.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "Poisson/PoissonDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<PureDiffFlux,
                       CellCenterFVMData,
                       ComputeDiffusiveFlux,
                       FiniteVolumePoissonModule>
PureDiffFluxProvider("PureDiffFlux");

//////////////////////////////////////////////////////////////////////////////

void PureDiffFlux::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

PureDiffFlux::PureDiffFlux(const std::string& name) :
  ComputeDiffusiveFlux(name),
  socket_volumes("volumes",false),
  _nbCVStates(0),
  _states(),
  _values(),
  _gradients(),
  _avState(),
  _varSet(CFNULL)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

PureDiffFlux::~PureDiffFlux()
{
  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void PureDiffFlux::computeFlux(RealVector& result)
{
  // reset the resulting flux to 0
  result = 0.0;

  GeometricEntity& geo = *getMethodData().getCurrentFace();
  SafePtr<DerivativeComputer> derivComputer = getMethodData().getDerivativeComputer();
  const bool isPerturb = this->getMethodData().isPerturb();
  
  if (!isPerturb) {
    // set the state values (pointers) corresponding to the vertices of the control volume
    derivComputer->computeControlVolume(_states, &geo);
    _nbCVStates = derivComputer->getNbVerticesInControlVolume(&geo);
    
  }
  _varSet->setGradientVars(_states, _values, _nbCVStates);
  
  // compute control volume around the face and gradients
  derivComputer->computeGradients(&geo, _values, _gradients);
  
  // compute the average values
  derivComputer->computeAverageValues(&geo, _states, _avState);
    
  const CFreal faceArea = this->socket_faceAreas.getDataHandle()[geo.getID()];

  // set the flux
  result = _varSet->getFlux(_avState, _gradients, getMethodData().getUnitNormal());
  result *= faceArea;

  CFLog(DEBUG_MED, "PureDiffFlux::computeFlux() => result = " << result << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void PureDiffFlux::setup()
{
  cf_assert(isConfigured());

  ComputeDiffusiveFlux::setup();

  _varSet = getMethodData().getDiffusiveVar().d_castTo<Physics::Poisson::PoissonDiffVarSet>();
  cf_assert(_varSet.isNotNull());

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  SafePtr<DerivativeComputer> derivComputer =
    getMethodData().getDerivativeComputer();
  
  // AL: careful here: this cannot assume that setup() has been run on derivComputer
  const CFuint nbNodesInControlVolume = derivComputer->getMaxNbVerticesInControlVolume();
  _states.resize(nbNodesInControlVolume);
  _values.resize(nbEqs, nbNodesInControlVolume);
  
  _gradients.resize(nbEqs);
  for (CFuint i = 0; i< nbEqs; ++i) {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }

  _avState.resize(PhysicalModelStack::getActive()->getNbEq());

  // set the diffusive flux jacobians to 0.
  _lFluxJacobian = 0.0;
  _rFluxJacobian = 0.0;
}
      
//////////////////////////////////////////////////////////////////////////////

void PureDiffFlux::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveFlux::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
PureDiffFlux::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
    ComputeDiffusiveFlux::needsSockets();
  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
