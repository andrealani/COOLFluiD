#include "LaxFriedCouplingFlux.hh"
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

MethodStrategyProvider<LaxFriedCouplingFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
laxFriedCouplingFluxProvider("LaxFriedCoupling");

//////////////////////////////////////////////////////////////////////////////

void LaxFriedCouplingFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >
    ("DiffCoeff","Diffusion reduction coefficient");
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedCouplingFlux::LaxFriedCouplingFlux(const std::string& name) :
  FVMCC_FluxSplitter(name),
  _sumFlux(),
  _aDiffVec(),
  _rightEv(),
  _leftEv(),
  _tempUnitNormal()
{
  addConfigOptionsTo(this);
  _currentDiffRedCoeff = 1.0;
  setParameter("DiffCoeff", &_currentDiffRedCoeff);
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedCouplingFlux::~LaxFriedCouplingFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedCouplingFlux::setup()
{
  FVMCC_FluxSplitter::setup();
  
  //  _avState = new State();
  _sumFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _aDiffVec.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());

}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedCouplingFlux::compute(RealVector& result)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  // flux for the right and left state
  _sumFlux = updateVarSet->getFlux()(pdata[1], unitNormal);
  _sumFlux += updateVarSet->getFlux()(pdata[0], unitNormal);
  
  updateVarSet->computeEigenValues(pdata[1], unitNormal, _rightEv);
  updateVarSet->computeEigenValues(pdata[0], unitNormal, _leftEv);
  
  (getMethodData().isPerturb()) ? computePerturbCase(result) : computeUnperturbCase(result);
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedCouplingFlux::computePerturbCase(RealVector& result)
{
  // here there could be the current number of equations if we want
  // different eigenvalues for each set of equations
  const EquationSubSysDescriptor& eqData =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  
  const CFuint nbEqs = eqData.getNbEqsSS();
  const CFuint start = eqData.getStartVarSS();
  
  CFreal aR = 0.0;
  CFreal aL = 0.0;
  CFuint count = start;
  for (CFuint i = 0; i < nbEqs; ++i, ++count) {
    aR = max(aR, std::abs(_rightEv[count]));
    aL = max(aL, std::abs(_leftEv[count]));
  }
  const CFreal a = max(aR,aL);

  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  vector<State*> *const solutionStates = getMethodData().getUpdateToSolutionVecTrans()->
    transformFromRefData(&pdata);
  
  // you must work with references (no copying allowed) !!!!
  State& leftState  = *(*solutionStates)[0];
  State& rightState = *(*solutionStates)[1];
  const CFreal aDiff = a*_currentDiffRedCoeff;

  result.slice(start, nbEqs) = 0.5*(_sumFlux.slice(start, nbEqs) - aDiff*
			     (rightState.slice(start, nbEqs) -
			      leftState.slice(start, nbEqs)));
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedCouplingFlux::computeUnperturbCase(RealVector& result)
{
  // here there could be the current number of equations if we want
  // different eigenvalues for each set of equations
  const EquationSubSysDescriptor& eqData =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
 
  const vector<vector<CFuint> >& eqPattern = eqData.getEqVarPatterns();
  const CFuint nbEqSS = eqData.getTotalNbEqSS();

  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  vector<State*> *const solutionStates = getMethodData().getUpdateToSolutionVecTrans()->
    transformFromRefData(&pdata);
  
  // you must work with references (no copying allowed) !!!!
  State& leftState  = *(*solutionStates)[0];
  State& rightState = *(*solutionStates)[1];
  
  for (CFuint iLSS = 0; iLSS < nbEqSS; ++iLSS) {
    const vector<CFuint>& eqIDs = eqPattern[iLSS];
    // unused // const CFuint start = eqIDs[0];
    const CFuint nbEqs = eqIDs.size();
    
    // detect the maximum eigenvalue in the current subset of equations
    CFreal aR = 0.0;
    CFreal aL = 0.0;
    for (CFuint i = 0; i < nbEqs; ++i) {
      aR = max(aR, std::abs(_rightEv[eqIDs[i]]));
      aL = max(aL, std::abs(_leftEv[eqIDs[i]]));
    }
    const CFreal aDiff  = max(aR,aL)*_currentDiffRedCoeff;
    for (CFuint i = 0; i < nbEqs; ++i) {
      _aDiffVec[eqIDs[i]] = aDiff;
    }
  }
  
  result = 0.5*(_sumFlux - _aDiffVec*(rightState - leftState));
  
  // compute update coefficient
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity& face = *data.getCurrentFace();
  const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
    data.getPolyReconstructor()->nbQPoints();
  
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  
  // compute eigen values for right state and opposite normal
  // only if the rigt state is not ghost
  if (!face.getState(1)->isGhost()) {
    _tempUnitNormal = -1.0*unitNormal;
    updateVarSet->computeEigenValues(pdata[1], _tempUnitNormal, _rightEv);
  }
  
  for (CFuint i = 0; i < nbEqSS; ++i) {
    // left contribution to update coefficient
    const CFuint leftID = face.getState(0)->getLocalID();
    const vector<CFuint>& eqIDs = eqPattern[i];
    const CFuint nbEqs = eqIDs.size();

    CFreal maxEigenVL = 0.0;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      maxEigenVL = max(_leftEv[eqIDs[iEq]], maxEigenVL);
    }
    updateCoeff[leftID*nbEqSS + i] += maxEigenVL*faceArea;

    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      CFreal maxEigenVR = 0.0;
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  maxEigenVR = max(_rightEv[eqIDs[iEq]], maxEigenVR);
	}

	const CFuint rightID = face.getState(1)->getLocalID();
	updateCoeff[rightID*nbEqSS + i] += maxEigenVR*faceArea;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
