#include "FiniteVolumeMHD/HLLDecEFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
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

MethodStrategyProvider<HLLDecEFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeMHDModule>
hllDecEFluxSplitterProvider("HLLDecE");

//////////////////////////////////////////////////////////////////////////////

void HLLDecEFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("UseAlpha", "Use diffusion reduction coefficient");
  options.addConfigOption< bool >("AddLax", "Add extra diffusion term to HLL solver");
  options.addConfigOption< CFreal >("Extra_kappaLax_E", "Amount of extra diffusion term added to HLL solver in Energy component");
  options.addConfigOption< bool >("2Dornot", "2D require to modify definition of plasma beta for extra diffusion");
}

//////////////////////////////////////////////////////////////////////////////

HLLDecEFlux::HLLDecEFlux(const std::string& name) :
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
  _addLax = false;
  setParameter("AddLax", &_addLax);
  _kappaLax_E = 0.125;
  setParameter("Extra_kappaLax_E", &_kappaLax_E);
  _2Dornot = false;
  setParameter("2Dornot", &_2Dornot);
}

//////////////////////////////////////////////////////////////////////////////

HLLDecEFlux::~HLLDecEFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void HLLDecEFlux::setup()
{
  FVMCC_FluxSplitter::setup();

  _rightFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void HLLDecEFlux::compute(RealVector& result)
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
     CFreal alpha=max(abs(amin),abs(amax))/(amax-amin);
	if (_addLax){
		alpha = 1.0;
	}
    result = 0.5*(_leftFlux + _rightFlux -
		  (amax + amin)/(amax - amin) * (_rightFlux - _leftFlux) +
		  alpha*2.0*amax * amin/(amax - amin) * (rightState - leftState));
  }

  //>> mark 2025.08 by HP
  //_addLax = true;
  if (_addLax){
	  GeometricEntity& face = *data.getCurrentFace();
	  State& left = *face.getState(0);
	  State& right = *face.getState(1);

	  CFreal P0 = 0.03851; // N/m**2
	  CFreal mu0 = 1.2566e-6;
	  CFreal p_tempL = left[7];
	  CFreal p_tempR = right[7];
	  CFreal pthL = p_tempL*P0;
	  CFreal pthR = p_tempR*P0;
	  CFreal BxL = left[4];
	  CFreal ByL = left[5];
	  CFreal BzL = left[6];
	  CFreal BxR = right[4];
	  CFreal ByR = right[5];
	  CFreal BzR = right[6];
	  CFreal Bsq2L = BxL*BxL + ByL*ByL + BzL*BzL;
	  CFreal Bsq2R = BxR*BxR + ByR*ByR + BzR*BzR;
	  CFreal pmagL = max(Bsq2L / 2. / mu0, 1.e-16);
	  CFreal pmagR = max(Bsq2R / 2. / mu0, 1.e-16);
	  CFreal plasmaBetaL = pthL / pmagL;
	  CFreal plasmaBetaR = pthR / pmagR;
	  CFreal fac = 1.0;
	  CFreal Q_BetafactorL = std::tanh(fac / plasmaBetaL / 100.0);
	  CFreal Q_BetafactorR = std::tanh(fac / plasmaBetaR / 100.0);
	  CFreal Q_Betafactor = max(Q_BetafactorL, Q_BetafactorR);
	  //>> Added BY Hp on 2025.10.14--------------
	  if (_2Dornot){
		  plasmaBetaL = p_tempL / max(Bsq2L / 2.0, 1.e-16);
		  plasmaBetaR = p_tempR / max(Bsq2R / 2.0, 1.e-16);
		  fac = 1.0;
		  Q_BetafactorL = std::tanh(fac / plasmaBetaL / 100.0);
		  Q_BetafactorR = std::tanh(fac / plasmaBetaR / 100.0);
		  Q_Betafactor = max(Q_BetafactorL, Q_BetafactorR);
		  Q_Betafactor = 1.0;
	  }
	  //<< ---------------------------------------
	  result[7] -= _kappaLax_E*Q_Betafactor*max(amax, std::abs(amin)) * (rightState[7] - leftState[7]);
  }
  //<< mark 2025.08 by HP
  
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
