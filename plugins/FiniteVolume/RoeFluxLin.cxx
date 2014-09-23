#include "RoeFluxLin.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/SubSystemStatus.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RoeFluxLin,
               BaseRoeFlux,
               FiniteVolumeModule,
               1>
linRoeFluxProvider("Linearized");

//////////////////////////////////////////////////////////////////////////////

RoeFluxLin::RoeFluxLin(const std::string& name) :
  BaseRoeFlux(name),
  _solutionVarSet(CFNULL),
  _updateVarSet(CFNULL),
  _linearizer(CFNULL),
  _solutionToLinearVarTrans(CFNULL),
  _updateToSolutionVarTrans(CFNULL),
  _sumFlux(),
  _rightEv(),
  _leftEv(),
  _eValues(),
  _rightEvalues(),
  _leftEvalues(),
  _absEvalues(),
  _jacobStar(),
  _absJacobStar(),
  _tempUnitNormal(),
  _solutionStates(CFNULL),
  _statesLR(2),
  _solutionStatesStar(CFNULL),
  _statesLRStar(2)
{
}

//////////////////////////////////////////////////////////////////////////////

RoeFluxLin::~RoeFluxLin()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& RoeFluxLin::compute()
{

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  FluxSplitterData& data = this->getFluxData();
  FluxSplitterData& dataStar = this->getFluxData2();

  dataStar.setPairID(data.getPairID());
  dataStar.setGeoShapeFunctions(data.getGeoShapeFunctions());

  _statesLR[0] = &data.getCurrLeftState();
  _statesLR[1] = &data.getCurrRightState();

  cf_assert(*_statesLR[0] == data.getCurrLeftState());
  cf_assert(*_statesLR[1] == data.getCurrRightState());

  _statesLRStar[0] = &dataStar.getCurrLeftState();
  _statesLRStar[1] = &dataStar.getCurrRightState();

  cf_assert(*_statesLRStar[0] == dataStar.getCurrLeftState());
  cf_assert(*_statesLRStar[1] == dataStar.getCurrRightState());

  _solutionStates = _updateToSolutionVarTrans->
    transformFromRef(&_statesLR);

  _solutionStatesStar = _updateToSolutionVarTrans->
    transformFromRef(&_statesLRStar);

  vector<State*> *const linearStatesStar = _solutionToLinearVarTrans->
    transform(_solutionStatesStar);

  _linearizer->linearize(*linearStatesStar);

//  std::cout << "data.unitNormal: " << data.unitNormal << std::endl;
  // set the eigenvectors and eigenvalues of the linearized jacobian
  _solutionVarSet->computeEigenValuesVectors(_rightEv,
                                         _leftEv,
                                         _eValues,
                                         data.unitNormal);

  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  setAbsEigenValues(dataStar);

  // abs of the jacobian
  _absJacobStar = _rightEv*(_absEvalues*_leftEv);
  _jacobStar = _rightEv*(_eValues*_leftEv);

  // flux for the right and left state
  State sumStatesStar = *_statesLRStar[0];
  sumStatesStar += *_statesLRStar[1];
  sumStatesStar *= 0.5;

  _sumFlux =  _solutionVarSet->getFlux()(sumStatesStar, dataStar.unitNormal);

  const State& stateL = *(*_solutionStates)[0];
  const State& stateR = *(*_solutionStates)[1];
  const State& stateLStar = *(*_solutionStatesStar)[0];
  const State& stateRStar = *(*_solutionStatesStar)[1];

  _result = _sumFlux - 0.5*(_jacobStar*(stateRStar + stateLStar)) ;
  _result += 0.5*(_jacobStar*(stateR + stateL)) - 0.5*(_absJacobStar*(stateR - stateL));

  // compute update coefficient
  if (!data.isPerturb) {

    DataHandle<CFreal> updateCoeff = data.updateCoeff->getDataHandle();
    const CFreal faceArea =
      (*data.faceAreas)[data.face->getID()]/data.nbQPoints;

    // left contribution to update coefficient
    CFreal maxEV =
      _updateVarSet->getMaxEigenValue(*_statesLR[0], data.unitNormal);

    const CFuint leftID = data.face->getState(0)->getLocalID();
    updateCoeff[leftID] += max(maxEV, 0.)*faceArea;

    if (!data.face->getState(1)->isGhost()) {
      // right contribution to update coefficient
      _tempUnitNormal = -1.0*data.unitNormal;
      maxEV = _updateVarSet->getMaxEigenValue(*_statesLR[1],_tempUnitNormal);
      const CFuint rightID = data.face->getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, 0.)*faceArea;
    }
  }

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& RoeFluxLin::computeJacobian()
{

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  FluxSplitterData& data = this->getFluxData();
  FluxSplitterData& dataStar = this->getFluxData2();

  dataStar.setPairID(data.getPairID());
  dataStar.setGeoShapeFunctions(data.getGeoShapeFunctions());

  _statesLR[0] = &data.getCurrLeftState();
  _statesLR[1] = &data.getCurrRightState();

  cf_assert(*_statesLR[0] == data.getCurrLeftState());
  cf_assert(*_statesLR[1] == data.getCurrRightState());

  _statesLRStar[0] = &dataStar.getCurrLeftState();
  _statesLRStar[1] = &dataStar.getCurrRightState();

  cf_assert(*_statesLRStar[0] == dataStar.getCurrLeftState());
  cf_assert(*_statesLRStar[1] == dataStar.getCurrRightState());

  _solutionStates = _updateToSolutionVarTrans->
    transformFromRef(&_statesLR);

  _solutionStatesStar = _updateToSolutionVarTrans->
    transformFromRef(&_statesLRStar);

  vector<State*> *const linearStatesStar = _solutionToLinearVarTrans->
    transform(_solutionStatesStar);

  _linearizer->linearize(*linearStatesStar);

//  std::cout << "data.unitNormal: " << data.unitNormal << std::endl;
  // set the eigenvectors and eigenvalues of the linearized jacobian
  _solutionVarSet->computeEigenValuesVectors(_rightEv,
                                         _leftEv,
                                         _eValues,
                                         data.unitNormal);

  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  setAbsEigenValues(dataStar);

  // abs of the jacobian
  _absJacobStar = _rightEv*(_absEvalues*_leftEv);
  _jacobStar = _rightEv*(_eValues*_leftEv);

  _matResult = 0.5*(_jacobStar) - 0.5*(_absJacobStar);

  return _matResult;
}

//////////////////////////////////////////////////////////////////////////////

void RoeFluxLin::setAbsEigenValues(FluxSplitterData& data)
{
  _absEvalues.abs(_eValues);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
