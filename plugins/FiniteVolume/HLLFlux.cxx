#include "HLLFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HLLFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
hllFluxSplitterProvider("HLL");

//////////////////////////////////////////////////////////////////////////////

void HLLFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("UseAlpha", "Use diffusion reduction coefficient");
}

//////////////////////////////////////////////////////////////////////////////

HLLFlux::HLLFlux(const std::string& name) :
  FVMCC_FluxSplitter(name),
  _rightFlux(),
  _leftFlux(),
  _rightEv(),
  _leftEv(),
  _tempUnitNormal(), 
  _solutionStates(CFNULL),
  _statesLR(2)
{
  addConfigOptionsTo(this);
  _useAlpha = false;
  setParameter("UseAlpha", &_useAlpha); 
}

//////////////////////////////////////////////////////////////////////////////

HLLFlux::~HLLFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void HLLFlux::setup()
{
  FVMCC_FluxSplitter::setup();

  _rightFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void HLLFlux::compute(RealVector& result)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  
  // flux for the right and left state
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  _rightFlux = updateVarSet->getFlux()(pdata[1], unitNormal);
  _leftFlux = updateVarSet->getFlux()(pdata[0], unitNormal);
  
  // here put update variables
  updateVarSet->computeEigenValues(pdata[1], unitNormal, _rightEv);
  updateVarSet->computeEigenValues(pdata[0], unitNormal, _leftEv);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  CFreal aR = 0.0;
  CFreal aL = 0.0;
  for (CFuint i = 0; i < nbEqs; ++i) {
    aR = max(aR, _rightEv[i]);
    aL = max(aL, _leftEv[i]);
  }

  const CFreal amax  = max((CFreal)0.0,max(aR,aL));

  aR = 0.0;
  aL = 0.0;
  for (CFuint i = 0; i < nbEqs; ++i) {
    aR = min(aR, _rightEv[i]);
    aL = min(aL, _leftEv[i]);
  }

  const CFreal amin  = min((CFreal)0.0,min(aR,aL));
  
  CellCenterFVMData& data = this->getMethodData(); 
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  _statesLR[0] = &polyRec->getCurrLeftState();
  _statesLR[1] = &polyRec->getCurrRightState();
  
  if (!getMethodData().reconstructSolVars()) {
    _solutionStates = getMethodData().getUpdateToSolutionVecTrans()->transform(&_statesLR);
  }
  else {
    _solutionStates = &_statesLR;
  }
  
  // you must work with references (no copying allowed) !!!!
  const State& leftState  = *(*_solutionStates)[0];
  const State& rightState = *(*_solutionStates)[1];

  if (!_useAlpha) {
    result = 0.5*(_leftFlux + _rightFlux -
		  (amax + amin)/(amax - amin) * (_rightFlux - _leftFlux) +
		  amax * amin/(amax - amin) * (rightState - leftState));
  }
  else {
    const CFreal alpha=max(abs(amin),abs(amax))/(amax-amin);
    result = 0.5*(_leftFlux + _rightFlux -
		  (amax + amin)/(amax - amin) * (_rightFlux - _leftFlux) +
		  alpha*2.0*amax * amin/(amax - amin) * (rightState - leftState));
  }
  
  // compute update coefficient
  if (!getMethodData().isPerturb()) {
    CellCenterFVMData& data = this->getMethodData(); 
    GeometricEntity& face = *data.getCurrentFace();
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      data.getPolyReconstructor()->nbQPoints();
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    // left contribution to update coefficient
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(_leftEv.max(), (CFreal)0.0)*faceArea;
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      const CFuint rightID = face.getState(1)->getLocalID();
      _tempUnitNormal = -1.0*unitNormal;
      const CFreal maxEV =
	updateVarSet->getMaxEigenValue(pdata[1], _tempUnitNormal);
      updateCoeff[rightID] += max(maxEV, (CFreal)0.0)*faceArea;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
