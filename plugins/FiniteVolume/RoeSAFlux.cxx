#include "RoeSAFlux.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"

#include <fstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RoeSAFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
roeSAFluxProvider("RoeSA");
      
//////////////////////////////////////////////////////////////////////////////

RoeSAFlux::RoeSAFlux(const std::string& name) :
  RoeFlux(name),
  _isPartitionFace(),
  _dataLeftState(),
  _dataRightState(),
  _tmpEv(),
  _normal()
{
  addConfigOptionsTo(this);
  _entropyFixID = 0;
  setParameter("entropyFixID",&_entropyFixID);
  
  _oldStencil = false;
  setParameter("oldStencil",&_oldStencil);
}

//////////////////////////////////////////////////////////////////////////////

RoeSAFlux::~RoeSAFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeSAFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("entropyFixID","ID of the entropy correction type.");
  options.addConfigOption< bool >
    ("oldStencil","Enlarged stencil used in old COOLFluiD, involving 2 more faces (8- instead of H-correction).");
}

//////////////////////////////////////////////////////////////////////////////

void RoeSAFlux::setAbsEigenValues()
{  
  //compute eigen values of the left state
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  const RealVector& unitNormal = getMethodData().getUnitNormal();
  
  updateVarSet->computeEigenValues(pdata[0], 
				   unitNormal,
				   _leftEvalues);
  
  //compute eigen values of the right state
  updateVarSet->computeEigenValues(pdata[1],
				   unitNormal,
				   _rightEvalues);
  
  //compute eigen values of the left state
  GeometricEntity& face = *getMethodData().getCurrentFace();
  CFreal etaSA = 0.0;

  // AL: the eigenvalue corresponding to the current face is computed twice: 
  // once from the left and once from the right => it is redundant 
  
  // 2D case
  // build the left cell
  updateEtaPA(face.getNeighborGeo(LEFT), _leftEvalues, etaSA);
  
  // build the right cell
  if (!face.getState(1)->isGhost()) {
    updateEtaPA(face.getNeighborGeo(RIGHT), _rightEvalues, etaSA);
  }
  etaSA *= 0.5;
  assert(std::abs(etaSA) > 0.);
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint i = 0; i < nbEqs; ++i) {
    const CFreal absEv = std::abs(_eValues[i]);
    
    switch (_entropyFixID) {
    case 1:
      _absEvalues[i] = (absEv >= 2.*etaSA) ? absEv : 0.25*absEv*absEv/etaSA + etaSA; 
      break;
    case 2:
      _absEvalues[i] = max(absEv, etaSA);
      break;
    case 3:
      _absEvalues[i] = absEv + etaSA;
      break;
    default:
      _absEvalues[i] = absEv + etaSA;
      break;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////
      
void RoeSAFlux::setup()
{
  RoeFlux::setup();
  
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_dataLeftState);
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(_dataRightState);
  _tmpEv.resize(PhysicalModelStack::getActive()->getNbEq());
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  
  getMethodData().setBuildAllCells(true);  
  
  Common::SafePtr<TopologicalRegionSet> pFaces = MeshDataStack::getActive()->
    getTrs("PartitionFaces");
  const CFuint nbPFaces = pFaces->getLocalNbGeoEnts();
  _isPartitionFace.resize(MeshDataStack::getActive()->Statistics().getNbFaces(), false);
  for (CFuint i = 0; i < nbPFaces; ++i) {
    _isPartitionFace[pFaces->getLocalGeoID(i)] = true;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void RoeSAFlux::updateEtaPA(GeometricEntity *const cell, 
			    const RealVector& eValues,
			    CFreal& etaSA) 
{
  GeometricEntity& currFace = *getMethodData().getCurrentFace();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  
  const GeomEntList *const faces = cell->getNeighborGeos();
  const CFuint nbFaces = faces->size();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    GeometricEntity *const neighborFace = (*faces)[iFace];
    // only faces sharing nodes with the orginal face are considered
    if (!_isPartitionFace[neighborFace->getID()]) {
      
      // the following check was not present in CF1 => the considered stencil is different (1 state less) 
      // => results are different
      if ((neighborFace->sharesNodesWith(currFace) && !_oldStencil) || _oldStencil) {      
	State *const rs = (neighborFace->getState(0) == cell->getState(0)) ? 
	  neighborFace->getState(1) : neighborFace->getState(0);
	assert(rs != NULL);
	
	const CFuint startID = neighborFace->getID()*dim;
	assert(neighborFace->getID() < socket_faceAreas.getDataHandle().size());
	const CFreal invArea = 1./socket_faceAreas.getDataHandle()[neighborFace->getID()];
	
	for (CFuint i = 0; i < dim; ++i) {
	  _normal[i] = normals[startID + i]*invArea;  
	}	
	
	updateVarSet->computePhysicalData(*rs, _dataRightState);  
	updateVarSet->computeEigenValues(_dataRightState,_normal, _tmpEv);
	for (CFuint i = 0; i < nbEqs; ++i) {
	  etaSA = max(etaSA, std::abs(_tmpEv[i] - eValues[i]));       
	}
      }
    }
  }
}	
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
