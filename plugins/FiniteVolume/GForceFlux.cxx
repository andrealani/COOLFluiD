#include "GForceFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GForceFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
stdGForceFluxProvider("GForce");

//////////////////////////////////////////////////////////////////////////////

GForceFlux::GForceFlux(const std::string& name) :
 FVMCC_FluxSplitter(name),
 _fluxR(),
 _fluxL(),
 _rightEv(),
 _leftEv(),
 _tempUnitNormal(),
 _gradient(),
 _bCoeffArray(),
 _stateLW(CFNULL),
 _statesLR(2),
 _pdataLW(),
 _coeff(1.0)
{
  addConfigOptionsTo(this);
  _maxCoeff = 1.0;
  setParameter("MaxCoeff", &_maxCoeff);

  _gradVarID = 0;
  setParameter("GradVarID", &_gradVarID);

  _freezeCoeffRes = -1.2;
  setParameter("FreezeCoeffRes", &_freezeCoeffRes);
}

//////////////////////////////////////////////////////////////////////////////

GForceFlux::~GForceFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void GForceFlux::setup()
{
  FVMCC_FluxSplitter::setup();

  _fluxR.resize(PhysicalModelStack::getActive()->getNbEq());
  _fluxL.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _leftEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());
  _gradient.resize(PhysicalModelStack::getActive()->getDim());
  _stateLW = new State();
  
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_pdataLW);
}

//////////////////////////////////////////////////////////////////////////////

void GForceFlux::unsetup()
{
  _fluxR.resize(0);
  _fluxL.resize(0);
  _rightEv.resize(0);
  _leftEv.resize(0);
  _tempUnitNormal.resize(0);
  _gradient.resize(0);
  _bCoeffArray.resize(0);
  deletePtr(_stateLW);
  _pdataLW.resize(0);
  
  FVMCC_FluxSplitter::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void GForceFlux::compute(RealVector& result)
{
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  
  // right and left fluxes
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  _fluxR = updateVarSet->getFlux()(pdata[1], unitNormal);
  _fluxL = updateVarSet->getFlux()(pdata[0], unitNormal);
  
  // all right and left eigenvalues
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

  // compute the coefficient for the blending if the flux
  // is not been perturbed
  const CFreal res = SubSystemStatusStack::getActive()->getResidual();

  // remove this from here
  if (_bCoeffArray.size() == 0) {
    _bCoeffArray.resize(socket_faceAreas.getDataHandle().size());
  }

  GeometricEntity& face = *getMethodData().getCurrentFace();
  
  if (!getMethodData().isPerturb() && res > _freezeCoeffRes) {
    //  State& currStateL = *data.face->getState(0);
    //     const CFreal lCoeff = computeBlendingCoeff(currStateL);

//     State& currStateR = *data.face->getState(1);
//     CFreal rCoeff = 0.0;
//     if (!currStateR.isGhost()) {
//       rCoeff = computeBlendingCoeff(currStateR);
//     }
//     else {
//       rCoeff = lCoeff;
//     }
//     // _coeff = lCoeff; //min(lCoeff, rCoeff);

//     _bCoeffArray[data.face->getID()] = min(lCoeff, rCoeff);
//     //lCoeff;
    
    _bCoeffArray[face.getID()] = _maxCoeff;

    // const CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
//     if (nbIter%200 == 0) {
//       std::string file = "coeff.dat." + StringOps::to_str(nbIter);
//       ofstream fout(file.c_str(), ios::app);
//       RealVector node = 0.5*(currStateR.getCoordinates() +
// 			     currStateL.getCoordinates());
//       fout << node << " " << _coeff << endl;
//     }
  }

  const CFreal om = _bCoeffArray[face.getID()]; // _coeff
  const CFreal oEminOm = 1. - om;
  
  vector<State*> *const solutionStates = getMethodData().getUpdateToSolutionVecTrans()->
    transformFromRefData(&pdata);
  
  // you must work with references (no copying allowed) !!!!
  State& leftState  = *(*solutionStates)[0];
  State& rightState = *(*solutionStates)[1];

  result = oEminOm*0.5*((_fluxR + _fluxL) - a*(rightState - leftState));

  const CFreal ovA = 1./a;
  *_stateLW = 0.5*((rightState + leftState) - ovA*(_fluxR - _fluxL));
  updateVarSet->computePhysicalData(*_stateLW, _pdataLW);
  
  result += om*updateVarSet->getFlux()(_pdataLW, unitNormal);
  
  // compute update coefficient
  if (!getMethodData().isPerturb()) {
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      getMethodData().getPolyReconstructor()->nbQPoints();
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    // left contribution to update coefficient
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(_leftEv.max(), (CFreal)0.0)*faceArea;
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      _tempUnitNormal = -1.0*unitNormal;
      const CFreal maxEV = updateVarSet->getMaxEigenValue(pdata[1], _tempUnitNormal);
      
      const CFuint rightID = face.getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, (CFreal)0.0)*faceArea;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void GForceFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
    ("MaxCoeff", "Max value for the blending coefficient");
  options.addConfigOption< CFuint >
    ("GradVarID", "ID of the variable use to compute gradient");
  options.addConfigOption< CFreal, Config::DynamicOption<> >
    ("FreezeCoeffRes",
     "Threshold residual controlling the freezing of the blending coeff");
}

//////////////////////////////////////////////////////////////////////////////

CFreal GForceFlux::computeBlendingCoeff(State& state)
{
  // FluxSplitterData& data = *getMethodData().getFluxSplitterData();

//   // create the current cell with pointer to the faces
//   SafePtr<GeometricEntityPool<CellTrsGeoBuilder> > cellBuilder =
//     getMethodData().getCellTrsGeoBuilder();
//   cellBuilder->getDataGE().idx = state.getLocalID();
//   GeometricEntity *const cell = cellBuilder->buildGE();

//   const GeomEntList *const faces = cell->getNeighborGeos();
//   const CFuint nbFaces = faces->size();

//   CFreal maxGrad = 0.0;
//   CFreal currGrad = 0.0;
//   for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
//     const GeometricEntity *const face = (*faces)[iFace];
//     State& otherState = (&state == face->getState(0)) ?
//       *face->getState(1) : *face->getState(0);

//     // const CFreal dr = MathFunctions::getDistance
//     //       (otherState.getCoordinates(), state.getCoordinates());
//     //     const CFreal grad =
//     //       std::abs((otherState[_gradVarID] - state[_gradVarID]))/dr;

//     const CFreal dw = otherState[_gradVarID] - state[_gradVarID];
//     _gradient = dw / (otherState.getCoordinates() - state.getCoordinates());
//     const CFreal grad =
//       std::abs(MathFunctions::innerProd(_gradient,data.unitNormal));

//     maxGrad = max(maxGrad, grad);
//     if (face->getID() == face.getID()) {
//       currGrad = grad;
//     }
//   }

//   maxGrad = (maxGrad > 0.0) ? maxGrad : MathTools::MathConsts::CFrealEps();
//   cf_assert(maxGrad > 0.0);

//   cellBuilder->releaseGE();

//   return min(_maxCoeff, 1./(1. + currGrad/maxGrad));
  return 0.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
