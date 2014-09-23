#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolumeNavierStokes/DiffFluxCarbuncleFixEuler.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DiffFluxCarbuncleFixEuler,
                       CellCenterFVMData,
		       ComputeDiffusiveFlux,
                       FiniteVolumeNavierStokesModule>
diffFluxCarbuncleFixEulerProvider("DiffFluxCarbuncleFixEuler");

//////////////////////////////////////////////////////////////////////////////

void DiffFluxCarbuncleFixEuler::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Eps","Epsilon for controlling dissipation.");
}

//////////////////////////////////////////////////////////////////////////////

DiffFluxCarbuncleFixEuler::DiffFluxCarbuncleFixEuler(const std::string& name) :
  ComputeDiffusiveFlux(name),
  socket_volumes("volumes",false),
  socket_fix_active("fix_active"),
  _nbCVStates(0),
  _states(),
  _values(),
  _gradients(),
  _avState(),
  _avPhysData(),
  _muS(0.)
{
  addConfigOptionsTo(this);
  
  _eps = 0.1;
  setParameter("Eps",&_eps);
}

//////////////////////////////////////////////////////////////////////////////

DiffFluxCarbuncleFixEuler::~DiffFluxCarbuncleFixEuler()
{
  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFluxCarbuncleFixEuler::computeFlux(RealVector& result)
{
  // reset the resulting flux to 0
  result = 0.0;
  
  GeometricEntity& geo = *getMethodData().getCurrentFace();
  SafePtr<DerivativeComputer> derivComputer = getMethodData().getDerivativeComputer();
  const bool isPerturb = this->getMethodData().isPerturb();
  DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
  
  if (!isPerturb) { 
    // set the state values (pointers) corresponding to the vertices of the control volume
    derivComputer->computeControlVolume(_states, &geo);
    _nbCVStates = derivComputer->getNbVerticesInControlVolume(&geo);
    
    static CFuint nbIter = 0;
    if (SubSystemStatusStack::getActive()->getNbIter() > nbIter) {
      nbIter = SubSystemStatusStack::getActive()->getNbIter();
      fix_active = 0.0;
    }
  }
  
  // set radients of velocity magnitude, temperature and Mach number    
  setGradientVars(_states, _values, _nbCVStates);
  
  // compute control volume around the face and gradients
  derivComputer->computeGradients(&geo, _values, _gradients);
  
  // compute the average values
  derivComputer->computeAverageValues(&geo, _states, _avState);
    
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  const RealVector& lData = pdata[LEFT];
  const RealVector& rData = pdata[RIGHT];
  
  const CFreal mL = lData[EulerTerm::V]/lData[EulerTerm::A];
  const CFreal mR = rData[EulerTerm::V]/rData[EulerTerm::A];
  
  const State& stateL = *geo.getState(0);
  const State& stateR = *geo.getState(1);
  
  getMethodData().getUpdateVar()->computePhysicalData(_avState, _avPhysData);
  
  // attempt
  _avPhysData[EulerTerm::VX] = 0.5*(lData[EulerTerm::VX] + rData[EulerTerm::VX]);
  _avPhysData[EulerTerm::VY] = 0.5*(lData[EulerTerm::VY] + rData[EulerTerm::VY]);
  _avPhysData[EulerTerm::V] = sqrt(_avPhysData[EulerTerm::VX]*_avPhysData[EulerTerm::VX] + _avPhysData[EulerTerm::VY]*_avPhysData[EulerTerm::VY]);
  _avPhysData[EulerTerm::A] = 0.5*(lData[EulerTerm::A] + rData[EulerTerm::A]);
  _avPhysData[EulerTerm::RHO] = 0.5*(lData[EulerTerm::RHO] + rData[EulerTerm::RHO]);
  
  bool isCrossedByShock = false;
  if ((mL >= 1. && mR < 1.) || (mL < 1. && mR >= 1.)) {
    CFreal dMachDcsi = 0; // gradient of Mach in the velocity direction
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      dMachDcsi += (*_gradients[2])[iDim]*_avPhysData[EulerTerm::VX+iDim];
    }
    if (dMachDcsi < 0.) {
      isCrossedByShock = true;
    }
  }
  
  // compute source term only if the cell is in the shock layer
  if (isCrossedByShock) {
    const CFreal uAvg = _avPhysData[EulerTerm::VX];
    const CFreal vAvg = _avPhysData[EulerTerm::VY];
    const CFreal velAvg = _avPhysData[EulerTerm::V];
    const CFreal rhoAvg = _avPhysData[EulerTerm::RHO];
    const CFreal cosD = uAvg/velAvg;
    const CFreal sinD = vAvg/velAvg;
    const RealVector& n = getMethodData().getUnitNormal();
    
    // only 2D for now
    // directional gradient of velocity in eta
    const CFreal dVdEta = -(*_gradients[0])[XX]*sinD + (*_gradients[0])[YY]*cosD;
    
    // eta component of the normal to the face
    const CFreal faceArea = this->socket_faceAreas.getDataHandle()[geo.getID()];
    
    // the muS is frozen during the iterative process
    if (!isPerturb) {
      // this choice can be not optimal
      // const CFreal lambdaMax = max(lData[EulerTerm::V] + lData[EulerTerm::A], rData[EulerTerm::V] + rData[EulerTerm::A]);
      // const CFreal lambdaMax = _avPhysData[EulerTerm::V] + _avPhysData[EulerTerm::A];
      const CFreal lambdaMax = std::abs(_avPhysData[EulerTerm::VX]*n[XX] + _avPhysData[EulerTerm::VY]*n[YY]) + _avPhysData[EulerTerm::A];
      
      // this choice can be not optimal
      const CFreal h = faceArea; //MathFunctions::getDistance(stateR.getCoordinates(), stateL.getCoordinates());
      _muS = _eps*h*lambdaMax*rhoAvg;
      fix_active[stateL.getLocalID()] = 1.;
      
      const CFreal diffUpdateCoeff = _muS/_avPhysData[EulerTerm::RHO]*faceArea*faceArea/(derivComputer->getControlVolume());
      DataHandle<CFreal> updateCoeff = this->socket_updateCoeff.getDataHandle();
      // left contribution to update coefficient
      const CFuint leftID = geo.getState(0)->getLocalID();
      updateCoeff[leftID] += diffUpdateCoeff;
      
      if (!geo.getState(1)->isGhost()) { 
	fix_active[stateR.getLocalID()] = 1.;
	// right contribution to update coefficient
	const CFuint rightID = geo.getState(1)->getLocalID();
	updateCoeff[rightID] += diffUpdateCoeff;
      }
    }
    
    const CFreal nEtaFace = -sinD*n[XX] + cosD*n[YY];
    const CFreal csiTerm = _muS*dVdEta*nEtaFace;
    result[1] = csiTerm*cosD*faceArea;
    result[2] = csiTerm*sinD*faceArea;
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffFluxCarbuncleFixEuler::setup()
{
  cf_assert(isConfigured());

  ComputeDiffusiveFlux::setup();
    
  const CFuint nbGradVars = 3;
  SafePtr<DerivativeComputer> derivComputer = getMethodData().getDerivativeComputer();
  
  // AL: careful here: this cannot assume that setup() has been run on derivComputer
  const CFuint nbNodesInControlVolume = derivComputer->getMaxNbVerticesInControlVolume();
  _states.resize(nbNodesInControlVolume);
  _values.resize(nbGradVars, nbNodesInControlVolume);
  
  // only gradients of velocity magnitude, temperature and Mach number are required 
  _gradients.resize(nbGradVars); 
  for (CFuint i = 0; i< nbGradVars; ++i) {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }
  
  _avState.resize(PhysicalModelStack::getActive()->getNbEq());
  getMethodData().getUpdateVar().d_castTo<EulerVarSet>()->getModel()->resizePhysicalData(_avPhysData);
  
  // set the diffusive flux jacobians to 0.
  _lFluxJacobian = 0.0;
  _rFluxJacobian = 0.0;
  
  const CFuint nbCells =
    MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();
  
  DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
  fix_active.resize(nbCells);
}
      
//////////////////////////////////////////////////////////////////////////////

void DiffFluxCarbuncleFixEuler::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveFlux::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffFluxCarbuncleFixEuler::setGradientVars(const vector<RealVector*>& states,
						RealMatrix& values,
						CFuint stateSize)
{  
  // We suppose to have u T variables here 
  // gradients of V, T, M
  const CFreal R = 287.046;
  const CFreal gamma = 1.4;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  for (CFuint i = 0; i < stateSize; ++i) {
    const RealVector& state = *states[i];
    CFreal V2 = 0.;
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      V2 += state[1+iDim]*state[1+iDim];
    }
    values(0,i) = sqrt(V2);                               // velocity magnitude
    values(1,i) = state[dim + 1];                        // temperature
    values(2,i) = values(0,i)/sqrt(gamma*R*values(1,i)); // Mach number 
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
DiffFluxCarbuncleFixEuler::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
    ComputeDiffusiveFlux::needsSockets();
  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
DiffFluxCarbuncleFixEuler::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result = 
    ComputeDiffusiveFlux::providesSockets();
  result.push_back(&socket_fix_active);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
