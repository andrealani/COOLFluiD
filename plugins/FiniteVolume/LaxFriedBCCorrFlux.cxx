#include "LaxFriedBCCorrFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LaxFriedBCCorrFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
stdLaxFriedBCCorrFluxProvider("LaxFriedBCCorr");

//////////////////////////////////////////////////////////////////////////////

LaxFriedBCCorrFlux::LaxFriedBCCorrFlux(const std::string& name) :
  FVMCC_FluxSplitter(name),
  _sumFlux(),
  _rightEv(),
  _leftEv(),
  _tempUnitNormal(),
  _avState(CFNULL),
  _avPdata(),
  _statesLR(2)
{
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedBCCorrFlux::~LaxFriedBCCorrFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedBCCorrFlux::setup()
{
  FVMCC_FluxSplitter::setup();

  _avState = new State();

  _sumFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());
  
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_avPdata);
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedBCCorrFlux::unsetup()
{
  deletePtr(_avState);

  FVMCC_FluxSplitter::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedBCCorrFlux::compute(RealVector& result)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  if (!getMethodData().getUseAverageFlux()) {
    // flux for the right and left state
    _sumFlux = updateVarSet->getFlux()(pdata[1],unitNormal);
    _sumFlux += updateVarSet->getFlux()(pdata[0],unitNormal);

    updateVarSet->computeEigenValues(pdata[1], unitNormal, _rightEv);
    updateVarSet->computeEigenValues(pdata[0], unitNormal, _leftEv);

    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

    CFreal aR = 0.0;
    CFreal aL = 0.0;
    for (CFuint i = 0; i < nbEqs; ++i) {
      aR = max(aR, std::abs(_rightEv[i]));
      aL = max(aL, std::abs(_leftEv[i]));
    }
    const CFreal a = max(aR,aL);

    vector<State*> *const solutionStates = getMethodData().getUpdateToSolutionVecTrans()->
      transformFromRefData(&pdata);
    
    // you must work with references (no copying allowed) !!!!
    const State& leftState  = *(*solutionStates)[0];
    const State& rightState = *(*solutionStates)[1];
    result = 0.5*(_sumFlux - a*(rightState - leftState));

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
  else {
    // average state
    (*_avState) = 0.5*(getMethodData().getPolyReconstructor()->getCurrLeftState() + 
		       getMethodData().getPolyReconstructor()->getCurrRightState());
    
    updateVarSet->computePhysicalData(*_avState, _avPdata);
    // flux for the right and left state
    _sumFlux = updateVarSet->getFlux()(_avPdata, unitNormal);
    
    // while computing the eigenvalues, use the left state the left state
    // and the average state as the right one (don't use the right one,
    // a ghost state, that could contain a negative temperature)
    
    // attempt
    updateVarSet->computeEigenValues(_avPdata, unitNormal, _rightEv);
    
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    CFreal aR = 0.0;
    CFreal aL = 0.0;
    for (CFuint i = 0; i < nbEqs; ++i) {
      aR = max(aR, std::abs(_rightEv[i]));
      aL = max(aL, std::abs(_leftEv[i]));
    }
    const CFreal a = max(aR,aL);

    vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
    vector<State*> *const solutionStates = getMethodData().getUpdateToSolutionVecTrans()->
      transformFromRefData(&pdata);

    // you must work with references (no copying allowed) !!!!
    const State& leftState  = *(*solutionStates)[0];
    const State& rightState = *(*solutionStates)[1];
    result = _sumFlux - 0.5*a*(rightState - leftState);

    // compute update coefficient
    if (!getMethodData().isPerturb()) {
      CellCenterFVMData& data = this->getMethodData(); 
      GeometricEntity& face = *data.getCurrentFace();
      const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
	data.getPolyReconstructor()->nbQPoints();
      
      // left contribution to update coefficient
      const CFuint leftID = face.getState(0)->getLocalID();
      
      // check this !!!!
      //  data.updateCoeff[leftID] += a*faceArea;
      // max(_leftEv.max(), 0.0)*faceArea;
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
      updateCoeff[leftID] += max(_leftEv.max(), (CFreal)0.0)*faceArea;
      
      // this treatment is only used on the boundary
      cf_assert(face.getState(1)->isGhost());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
