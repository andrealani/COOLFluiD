#include "RoeFluxALE.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
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

MethodStrategyProvider<RoeFluxALE,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
roeFluxALEProvider("RoeALE");
      
//////////////////////////////////////////////////////////////////////////////

void RoeFluxALE::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("entropyFixID","ID of the entropy correction type.");
}

//////////////////////////////////////////////////////////////////////////////

RoeFluxALE::RoeFluxALE(const std::string& name) :
  RoeFlux(name),
  socket_pastNodes("pastNodes"),
  socket_futureNodes("futureNodes"),
  _vgn(0.),
  _meshSpeed()  
{
  addConfigOptionsTo(this);
  
  _entropyFixID = 0;
  setParameter("entropyFixID",&_entropyFixID);
}

//////////////////////////////////////////////////////////////////////////////

RoeFluxALE::~RoeFluxALE()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeFluxALE::compute(RealVector& result)
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
  
  if (!getMethodData().reconstructSolVars()) {
    _solutionStates = getMethodData().getUpdateToSolutionVecTrans()->transform(&_statesLR);
  }
  else {
    _solutionStates = &_statesLR;
  }
  
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
  
  // this only works with conservative variables
  _sumFlux -= (*_statesLR[0]) * _vgn;
  _sumFlux -= (*_statesLR[1]) * _vgn;
  
  result = 0.5*(_sumFlux - _absJacob*(stateR - stateL));

  // compute update coefficient
 if (!getMethodData().isPerturb()) {
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      polyRec->nbQPoints();
    
    // left contribution to update coefficient
    CFreal maxEV = updateVarSet->getMaxEigenValue(pdata[0], unitNormal);
    
    //Modify the maximum eigen values (due to mesh movement)
    maxEV -= _vgn;
    
    // left contribution to update coefficient
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(maxEV, (CFreal)0.)*faceArea;
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      
      _tempUnitNormal = -1.0*unitNormal;
      maxEV = updateVarSet->getMaxEigenValue(pdata[1],_tempUnitNormal);
      
      //Modify the maximum eigen values (due to mesh movement)
      maxEV -= (-1.0)*_vgn;
      
      const CFuint rightID = face.getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, (CFreal)0.)*faceArea;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RoeFluxALE::setAbsEigenValues()
{
  //Modify the eigen values to account for mesh deformation
  /// Compute the meshSpeed
  CFreal dt = SubSystemStatusStack::getActive()->getDT();
  
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity *const geo = data.getCurrentFace();
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  const RealVector& shapeFunction = polyRec->getCurrentGeoShapeFunction(geo);
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  
  _meshSpeed = 0.;
  
  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();
  
  const CFuint nbNodes = geo->nbNodes();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  // Compute speed of the mesh at current quadrature point
  for(CFuint iNode = 0; iNode < nbNodes; iNode++){
    const CFuint nodeID = (geo->getNode(iNode))->getLocalID();
    for(CFuint iDim = 0; iDim < nbDim; iDim++){
      // nodal shape function value in correspondance of the quadrature point position
      _meshSpeed[iDim] += shapeFunction[iNode] * (*(futureNodes[nodeID]))[iDim];
      _meshSpeed[iDim] -= shapeFunction[iNode] * (*(pastNodes[nodeID]))[iDim];
    }
  }

  _meshSpeed /= dt;
  
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  //Compute vg*n
  _vgn = _meshSpeed[0] * unitNormal[0];
  for(CFuint iDim = 1;iDim < nbDim ;iDim++){
    _vgn += _meshSpeed[iDim] * unitNormal[iDim];
  }
   
  vector<RealVector>& pdata = polyRec->getExtrapolatedPhysicaData();

  //compute eigen values of the left state
  updateVarSet->computeEigenValues(pdata[0], unitNormal, _leftEvalues);
  
  //compute eigen values of the right state
  updateVarSet->computeEigenValues(pdata[1], unitNormal, _rightEvalues);
  
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
      _absEvalues[i] = (!(std::abs(_eValues[i]) < 2.0)) ?
	absLambda : (absLambda*absLambda/(4.0*lambdaCorr) + lambdaCorr);
    }
  }
  if (_entropyFixID == 0) {
    //Set absolute value of eigenvalues
    _absEvalues = abs(_eValues);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > RoeFluxALE::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = RoeFlux::needsSockets();
  
  result.push_back(&socket_pastNodes);
  result.push_back(&socket_futureNodes);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RoeFluxALE::setup()
{
  RoeFlux::setup();
  _meshSpeed.resize(Framework::PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
