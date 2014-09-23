#include "RoeFluxLinALE.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
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

Environment::ObjectProvider<RoeFluxLinALE,
               BaseRoeFlux,
               FiniteVolumeModule,
               1>
RoeFluxLinALEProvider("LinearizedALE");

//////////////////////////////////////////////////////////////////////////////

void RoeFluxLinALE::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("entropyFixID","ID of the entropy correction type.");
}

//////////////////////////////////////////////////////////////////////////////

RoeFluxLinALE::RoeFluxLinALE(const std::string& name) :
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
  _meshSpeed(),
  _solutionStates(CFNULL),
  _statesLR(2),
  _solutionStatesStar(CFNULL),
  _statesLRStar(2),
  _pastNodes(CFNULL),
  _futureNodes(CFNULL)
{

  addConfigOptionsTo(this);

  _entropyFixID = 0;
  setParameter("entropyFixID",&_entropyFixID);

}

//////////////////////////////////////////////////////////////////////////////

RoeFluxLinALE::~RoeFluxLinALE()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& RoeFluxLinALE::compute()
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
    const CFreal faceArea =
      (*data.faceAreas)[data.face->getID()]/data.nbQPoints;

    // left contribution to update coefficient
    CFreal maxEV =
      _updateVarSet->getMaxEigenValue(*_statesLR[0], data.unitNormal);

    //Modify the maximum eigen values (due to mesh movement)
    maxEV -= _vgn;

    DataHandle<CFreal> updateCoeff = data.updateCoeff->getDataHandle();
    // left contribution to update coefficient
    const CFuint leftID = data.face->getState(0)->getLocalID();
    updateCoeff[leftID] += max(maxEV, 0.)*faceArea;

    if (!data.face->getState(1)->isGhost()) {
      // right contribution to update coefficient
      _tempUnitNormal = -1.0*data.unitNormal;
      maxEV = _updateVarSet->getMaxEigenValue(*_statesLR[1],_tempUnitNormal);

      //Modify the maximum eigen values (due to mesh movement)
      maxEV -= (-1.0)*_vgn;

      const CFuint rightID = data.face->getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, 0.)*faceArea;
    }
  }

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& RoeFluxLinALE::computeJacobian()
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

void RoeFluxLinALE::setAbsEigenValues(FluxSplitterData& data)
{
  //Modify the eigen values to account for mesh deformation
  /// Compute the meshSpeed
  CFreal dt = SubSystemStatusStack::getActive()->getDT();

  const RealVector shapeFunction = data.getGeoShapeFunction();
  GeometricEntity *const geo = data.face;

  _meshSpeed = 0.;

  const CFuint nbNodes = geo->nbNodes();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  // Compute speed of the mesh at current quadrature point
  for(CFuint iNode = 0; iNode < nbNodes; iNode++){

    const CFuint nodeID = (geo->getNode(iNode))->getLocalID();
    for(CFuint iDim = 0; iDim < nbDim; iDim++){
      _meshSpeed[iDim] += shapeFunction[iNode] * (*(_futureNodes[nodeID]))[iDim];
      _meshSpeed[iDim] -= shapeFunction[iNode] * (*(_pastNodes[nodeID]))[iDim];
    }
  }

  _meshSpeed /= dt;

  //Compute vg*n
  _vgn = _meshSpeed[0] * data.unitNormal[0];
  for(CFuint iDim = 1;iDim < nbDim ;iDim++){
    _vgn += _meshSpeed[iDim] * data.unitNormal[iDim];
  }

  //compute eigen values of the left state
  State& stateL = data.getCurrLeftState();
  State& stateR = data.getCurrRightState();

  _solutionVarSet->computeEigenValues(stateL,
				  data.unitNormal,
				  _leftEvalues);

  //compute eigen values of the right state
  _solutionVarSet->computeEigenValues(stateR,
				  data.unitNormal,
				  _rightEvalues);

  //Modify the eigen values (eValue = eValue - Vg.n)
  _eValues -= _vgn;
  _rightEvalues -= _vgn;
  _leftEvalues -= _vgn;

  CFreal lambdaCorr = 0.0;
  computeLambdaCorr(lambdaCorr);

  lambdaCorr *= 0.5;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  if (_entropyFixID == 1) {
    for (CFuint i = 0; i < nbEqs; ++i) {
      _absEvalues[i] = max(std::abs(_eValues[i]), lambdaCorr);
    }
  }
  else if (_entropyFixID == 2) {
    for (CFuint i = 0; i < nbEqs; ++i) {
      const CFreal absLambda = std::abs(_eValues[i]);
      _absEvalues[i] = (!std::abs(_eValues[i]) < 2.0) ?
	absLambda : (absLambda*absLambda/(4.0*lambdaCorr) + lambdaCorr);
    }
  }
  if (_entropyFixID == 0) {
    //Set absolute value of eigenvalues
    _absEvalues.abs(_eValues);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
