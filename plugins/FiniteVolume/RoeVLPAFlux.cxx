#include "RoeVLPAFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/BaseTerm.hh"
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

MethodStrategyProvider<RoeVLPAFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
roeVLPAFluxProvider("RoeVLPA");

//////////////////////////////////////////////////////////////////////////////

RoeVLPAFlux::RoeVLPAFlux(const std::string& name) :
  RoeFlux(name),
  _dataLeftState(),
  _dataRightState(),
  _tmpEv(),
  _normal()
{
  addConfigOptionsTo(this);
  _useFixPA = false;
  setParameter("useFixPA",&_useFixPA);

  _vlFixIDs = vector<CFuint>();
  setParameter("VLFixIDs",&_vlFixIDs);
  
  _flagVL = vector<bool>();
}
      
//////////////////////////////////////////////////////////////////////////////

RoeVLPAFlux::~RoeVLPAFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeVLPAFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("useFixPA","Flag to apply Pandolfi-D'Ambrosio's fix");
  options.addConfigOption< vector<CFuint> >
    ("VLFixIDs", "IDs of the variables on which VL fix must be applied");
}

//////////////////////////////////////////////////////////////////////////////

void RoeVLPAFlux::setAbsEigenValues()
{
  //compute eigen values of the left state
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  updateVarSet->computeEigenValues(pdata[0],
				   unitNormal,
				   _leftEvalues);
  
  //compute eigen values of the right state
  updateVarSet->computeEigenValues(pdata[1],
				   unitNormal,
				   _rightEvalues);
  
  
  CFreal etaVL = 0.0;
  computeEtaVL(etaVL);
  
  //compute eigen values of the left state
  GeometricEntity& face = *getMethodData().getCurrentFace();
  State *const lState = face.getState(LEFT);
  State *const rState = face.getState(RIGHT);
  CFreal etaPA = 0.0;
  
  updateEtaPA(face.getNeighborGeo(LEFT), _leftEvalues, rState, etaPA);
  
  // build the right cell
  if (!rState->isGhost()) {
    updateEtaPA(face.getNeighborGeo(RIGHT), _rightEvalues, lState, etaPA);
  }
  // is this needed ??
  etaPA *= 0.5;
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint i = 0; i < nbEqs; ++i) {
    const CFreal absEv = std::abs(_eValues[i]);
    
    if (_flagVL[i]) {
      if (!_useFixPA) {
	_absEvalues[i] = (absEv >= 2.*etaVL) ? absEv : 0.25*absEv*absEv/etaVL + etaVL;
      }
      else { 
	_absEvalues[i] = absEv;
      }
    }
    else {
      _absEvalues[i] = max(absEv, etaPA);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RoeVLPAFlux::setup()
{
  RoeFlux::setup();
  
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_dataLeftState);
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_dataRightState);
 
  _tmpEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  
  _flagVL.resize(PhysicalModelStack::getActive()->getNbEq());
  _flagVL.assign(_flagVL.size(), false);
  for (CFuint i = 0; i < _vlFixIDs.size(); ++i) {
    _flagVL[_vlFixIDs[i]] = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void RoeVLPAFlux::updateEtaPA(GeometricEntity *const cell, 
			      const RealVector& eValues,
			      State *const otherState,
			      CFreal& etaPA) 
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  
  const GeomEntList *const faces = cell->getNeighborGeos();
  const CFuint nbFaces = faces->size();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    const GeometricEntity *const neighborFace = (*faces)[iFace];
    State *const rs = (neighborFace->getState(0) == cell->getState(0)) ? 
      neighborFace->getState(1) : neighborFace->getState(0);
    
    //    if (rs != otherState || (rs == otherState && otherState->isGhost()) ) {
    if (rs != otherState) {
      const CFuint startID = neighborFace->getID()*dim;
      const CFreal invArea = 1./socket_faceAreas.getDataHandle()[neighborFace->getID()];
      for (CFuint i = 0; i < dim; ++i) {
	_normal[i] = normals[startID + i]*invArea;  
      }	      
      
      updateVarSet->computePhysicalData(*rs, _dataRightState);
      updateVarSet->computeEigenValues(_dataRightState,_normal, _tmpEv);
      
      for (CFuint i = 0; i < nbEqs; ++i) {
	etaPA = max(etaPA, std::abs(_tmpEv[i] - eValues[i]));       
      }
    }
  }
}	
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
