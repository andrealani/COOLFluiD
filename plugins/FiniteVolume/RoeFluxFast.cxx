#include "RoeFluxFast.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RoeFluxFast,
               BaseRoeFlux,
               FiniteVolumeModule,
               1>
fastRoeFluxProvider("Fast");

//////////////////////////////////////////////////////////////////////////////

RoeFluxFast::RoeFluxFast(const std::string& name) :
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
  _absJacob(),
  _tempUnitNormal(),
  _solutionStates(CFNULL),
  _statesLR(2)
{
}

//////////////////////////////////////////////////////////////////////////////

RoeFluxFast::~RoeFluxFast()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& RoeFluxFast::compute()
{
  FluxSplitterData& data = this->getFluxData();

  _statesLR[0] = &data.getCurrLeftState();
  _statesLR[1] = &data.getCurrRightState();

  cf_assert(*_statesLR[0] == data.getCurrLeftState());
  cf_assert(*_statesLR[1] == data.getCurrRightState());

  _solutionStates = _updateToSolutionVarTrans->
    transformFromRef(&_statesLR);

  vector<State*> *const linearStates = _solutionToLinearVarTrans->
    transform(_solutionStates);

  _linearizer->linearize(*linearStates);

  // set the eigenvectors and eigenvalues of the linearized jacobian
  _solutionVarSet->computeEigenValuesVectors(_rightEv,
                                         _leftEv,
                                         _eValues,
                                         data.unitNormal);

  // set the abs of the  eigen values (the implementation of this
  // function change if there are entropy or carbuncle fixes)
  setAbsEigenValues(data);
  // abs of the jacobian
  _absJacob = _rightEv*(_absEvalues*_leftEv);

  // flux for the right and left state
  State& sumStates = *_statesLR[0];
  sumStates += *_statesLR[1];
  sumStates *= 0.5;
  _sumFlux =  _solutionVarSet->getFlux()(sumStates, data.unitNormal);

  const State& stateL = *(*_solutionStates)[0];
  const State& stateR = *(*_solutionStates)[1];
  _result = _sumFlux - 0.5*(_absJacob*(stateR - stateL));

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

void RoeFluxFast::setAbsEigenValues(FluxSplitterData& data)
{
  _absEvalues.abs(_eValues);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
