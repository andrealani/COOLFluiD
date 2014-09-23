#include "RoeFluxTurb.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RoeFluxTurb,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
roeFluxTurbProvider("RoeTurb");

//////////////////////////////////////////////////////////////////////////////

RoeFluxTurb::RoeFluxTurb(const std::string& name) :
  RoeFlux(name)
{
}

//////////////////////////////////////////////////////////////////////////////

RoeFluxTurb::~RoeFluxTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeFluxTurb::compute(RealVector& result)
{
  (this->getMethodData().isPerturb()) ? 
    computePerturbCase(result) : computeUnperturbCase(result);
}
      
//////////////////////////////////////////////////////////////////////////////

void RoeFluxTurb::computePerturbCase(RealVector& result)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  SafePtr<ConvectiveVarSet> solutionVarSet = getMethodData().getSolutionVar();
  CellCenterFVMData& data = this->getMethodData(); 
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  
  _statesLR[0] = &polyRec->getCurrLeftState();
  _statesLR[1] = &polyRec->getCurrRightState();
  
  cf_assert(*_statesLR[0] == polyRec->getCurrLeftState());
  cf_assert(*_statesLR[1] == polyRec->getCurrRightState());
  
  _solutionStates = getMethodData().getUpdateToSolutionVecTrans()->transform(&_statesLR);
  
  linearize();
  
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  // set the eigenvectors and eigenvalues of the linearized jacobian
  solutionVarSet->computeEigenValuesVectors(_rightEv,
					    _leftEv,
					    _eValues,
					    unitNormal);
  
  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  setAbsEigenValues();
  
  // abs of the jacobian
  _absJacob = _rightEv*(_absEvalues*_leftEv);

  // flux for the right and left state
  vector<RealVector>& pdata = polyRec->getExtrapolatedPhysicaData();
  _sumFlux =  updateVarSet->getFlux()(pdata[0], unitNormal);
  _sumFlux += updateVarSet->getFlux()(pdata[1], unitNormal);
  
  const State& stateL = *(*_solutionStates)[0];
  const State& stateR = *(*_solutionStates)[1];

  const EquationSubSysDescriptor& eqData =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();
  //  const CFuint start = eqData.getStartVarSS();
  
  //   const vector<vector<CFuint> >& eqPattern = eqData.getEqVarPatterns();
  //   const CFuint nbEqSS = eqData.getTotalNbEqSS();
  
  ///Roe for every variable
   result = 0.5*(_sumFlux - _absJacob*(stateR - stateL));
   
   ///LaxFriedrichs for Turbulent quantities
   solutionVarSet->computeEigenValues(pdata[1], unitNormal, _rightEvalues);
   solutionVarSet->computeEigenValues(pdata[0], unitNormal, _leftEvalues);
   
   CFreal aR = 0.0;
   CFreal aL = 0.0;
   for (CFuint i = 4; i < nbEqs; ++i) {
     aR = max(aR, std::abs(_rightEvalues[i]));
     aL = max(aL, std::abs(_leftEvalues[i]));
   }
   const CFreal a = max(aR,aL);
   const CFreal aDiff = a*_currentDiffRedCoeff;
   
   for (CFuint iEq = 4; iEq < nbEqs; ++iEq) {
     result[iEq] = 0.5*(_sumFlux[iEq] - aDiff*(stateR[iEq] - stateL[iEq]));
   }

// if((nbEqs >= 4) && (start == 0)){
//   result[4] = _unperturbedFlux[4];
//   result[5] = _unperturbedFlux[5];
// }
// else{
//   result[0] = _unperturbedFlux[0];
//   result[1] = _unperturbedFlux[1];
//   result[2] = _unperturbedFlux[2];
//   result[3] = _unperturbedFlux[3];
//
//   CFuint currentVar = data.iPerturbVar;
//   //if 2 turb equations coupled
//   if(nbEqs == 1){
//     if(currentVar == 4) result[5] = _unperturbedFlux[5];
//     if(currentVar == 5) result[4] = _unperturbedFlux[4];
//   }
//   if(nbEqs == 2){
//   }
// }

//   FluxSplitterData& data = this->getFluxData();
//
//   const EquationSubSysDescriptor& eqData =
//     PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
//   const CFuint nbEqs = eqData.getNbEqsSS();
//   const CFuint start = eqData.getStartVarSS();
//
//   _statesLR[0] = &data.getCurrLeftState();
//   _statesLR[1] = &data.getCurrRightState();
//
//   cf_assert(*_statesLR[0] == data.getCurrLeftState());
//   cf_assert(*_statesLR[1] == data.getCurrRightState());
//
//   _solutionStates = _updateToSolutionVarTrans->
//     transformFromRef(&_statesLR);
//
//   // flux for the right and left state
//   _sumFlux =  _solutionVarSet->getFlux()(*_statesLR[0], data.unitNormal);
//   _sumFlux += _solutionVarSet->getFlux()(*_statesLR[1], data.unitNormal);
//
//   State& stateL = *(*_solutionStates)[0];
//   State& stateR = *(*_solutionStates)[1];
//
//   RealSliceVector::setSize(nbEqs);
//   RealSliceMatrix::setNbRowsCols(nbEqs,nbEqs);
//
// if((nbEqs >= 4) && (start == 0)){
//
//   vector<State*> *const linearStates = _solutionToLinearVarTrans->
//     transform(_solutionStates);
//
//   _linearizer->linearize(*linearStates);
//
//   // set the eigenvectors and eigenvalues of the linearized jacobian
//   _solutionVarSet->computeEigenValuesVectors(_rightEv,
//                                          _leftEv,
//                                          _eValues,
//                                          data.unitNormal);
//
//   // set the abs of the  eigen values (the implementation of this
//   // function change if there are entropy or carbuncle fixes)
//   setAbsEigenValues(data);
//
//   // abs of the jacobian
//   _absJacob.slice(start,start) =   _rightEv.slice(start,start)*(_absEvalues.slice(start)*_leftEv.slice(start,start));
//
//   ///Roe for every variable
//   result.slice(start) =
// 0.5*(_sumFlux.slice(start) - _absJacob.slice(start,start)*(stateR.slice(start) - stateL.slice(start)));
//
// }
// else{
// //if 2 turb equations coupled
// if(nbEqs == 2){
//   ///LaxFriedrichs for Turbulent quantities
//
//   _solutionVarSet->computeEigenValues(data.getCurrRightState(),
//    				  data.unitNormal, _rightEvalues);
//   _solutionVarSet->computeEigenValues(data.getCurrLeftState(),
//    				  data.unitNormal, _leftEvalues);
//
//   CFreal aR = 0.0;
//   CFreal aL = 0.0;
//   for (CFuint i = 0; i < nbEqs; ++i) {
//     aR = max(aR, std::abs(_rightEvalues[start+i]));
//     aL = max(aL, std::abs(_leftEvalues[start+i]));
//   }
//   const CFreal a = max(aR,aL);
//   const CFreal aDiff = a*_currentDiffRedCoeff;
//
//   result.slice(start) = 0.5*(_sumFlux.slice(start) - aDiff*(stateR.slice(start) - stateL.slice(start)));
//
// }
// //if 2 turb equations are uncoupled
// if(nbEqs == 1){
//   ///LaxFriedrichs for Turbulent quantities
//   _solutionVarSet->computeEigenValues(data.getCurrRightState(),
//    				  data.unitNormal, _rightEvalues);
//   _solutionVarSet->computeEigenValues(data.getCurrLeftState(),
//    				  data.unitNormal, _leftEvalues);
//
//   CFreal aR = 0.0;
//   CFreal aL = 0.0;
//   for (CFuint i = 0; i < nbEqs; ++i) {
//     aR = max(aR, std::abs(_rightEvalues[start+i]));
//     aL = max(aL, std::abs(_leftEvalues[start+i]));
//   }
//   const CFreal a = max(aR,aL);
//   const CFreal aDiff = a*_currentDiffRedCoeff;
//
//   result.slice(start) = 0.5*(_sumFlux.slice(start) - aDiff*(stateR.slice(start) - stateL.slice(start)));
// }
//
// }

}

//////////////////////////////////////////////////////////////////////////////

void RoeFluxTurb::computeUnperturbCase(RealVector& result)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  SafePtr<ConvectiveVarSet> solutionVarSet = getMethodData().getSolutionVar();
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity& face = *data.getCurrentFace();
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  
  _statesLR[0] = &polyRec->getCurrLeftState();
  _statesLR[1] = &polyRec->getCurrRightState();
  
  cf_assert(*_statesLR[0] == polyRec->getCurrLeftState());
  cf_assert(*_statesLR[1] == polyRec->getCurrRightState());
  
  _solutionStates = getMethodData().getUpdateToSolutionVecTrans()->transform(&_statesLR);
  
  linearize();
  
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  // set the eigenvectors and eigenvalues of the linearized jacobian
  solutionVarSet->computeEigenValuesVectors(_rightEv,
					    _leftEv,
					    _eValues,
					    unitNormal);
  
  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  setAbsEigenValues();
  // abs of the jacobian
  _absJacob = _rightEv*(_absEvalues*_leftEv);
  
   // flux for the right and left state
  vector<RealVector>& pdata = polyRec->getExtrapolatedPhysicaData();
  _sumFlux =  updateVarSet->getFlux()(pdata[0], unitNormal);
  _sumFlux += updateVarSet->getFlux()(pdata[1], unitNormal);
  
  const State& stateL = *(*_solutionStates)[0];
  const State& stateR = *(*_solutionStates)[1];
  
  const EquationSubSysDescriptor& eqData =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();
  //unused//  const CFuint start = eqData.getStartVarSS();
  const vector<vector<CFuint> >& eqPattern = eqData.getEqVarPatterns();
  const CFuint nbEqSS = eqData.getTotalNbEqSS();
  
  // Roe for every variable
  result = 0.5*(_sumFlux - _absJacob*(stateR - stateL));
  
  //LaxFriedrichs for Turbulent quantities
  solutionVarSet->computeEigenValues(pdata[1], unitNormal, _rightEvalues);
  solutionVarSet->computeEigenValues(pdata[0], unitNormal, _leftEvalues);
  
  CFreal aR = 0.0;
  CFreal aL = 0.0;
  for (CFuint i = 4; i < nbEqs; ++i) {
    aR = max(aR, std::abs(_rightEvalues[i]));
    aL = max(aL, std::abs(_leftEvalues[i]));
  }
  const CFreal a = max(aR,aL);
  const CFreal aDiff = a*_currentDiffRedCoeff;

  for (CFuint iEq = 4; iEq < nbEqs; ++iEq) {
    result[iEq] = 0.5*(_sumFlux[iEq] - aDiff*(stateR[iEq] - stateL[iEq]));
  }
  
  // compute update coefficient
  const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
     polyRec->nbQPoints();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // compute eigen values for right state and opposite normal
  // only if the rigt state is not ghost
  if (!face.getState(1)->isGhost()) {
    _tempUnitNormal = -1.0*unitNormal;
    solutionVarSet->computeEigenValues(pdata[1], _tempUnitNormal, _rightEvalues);
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
	maxEigenVR = max(_rightEvalues[eqIDs[iEq]], maxEigenVR);
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
