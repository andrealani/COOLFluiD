#include "FiniteVolume/FiniteVolume.hh"
#include "ShiftedPeriodicX2D.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ShiftedPeriodicX2D, CellCenterFVMData, FiniteVolumeModule> shiftedPeriodicX2DFVMCCProvider("ShiftedPeriodicX2DFVMCC");
    
//////////////////////////////////////////////////////////////////////////////

ShiftedPeriodicX2D::ShiftedPeriodicX2D(const std::string& name) : 
  FVMCC_BC(name),
  _periodicFaceID(),
  _boundaryStatesInSequence(),
  _localTRSFaceIDInSequence()	
{
}

//////////////////////////////////////////////////////////////////////////////

ShiftedPeriodicX2D::~ShiftedPeriodicX2D() 
{
}

//////////////////////////////////////////////////////////////////////////////

void ShiftedPeriodicX2D::setup()
{
  FVMCC_BC::setup();
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  
  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
	      geoBuilder = getMethodData().getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  geoData.trs = trs;

  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();

  _periodicFaceID.resize(nbTrsFaces);
  _boundaryStatesInSequence.resize(nbTrsFaces);

  // mapping the global face IDs to local IDs in the TRS and
  // storing the corresponding boundary states in a temporary vector
  const CFuint nbTrFaces = nbTrsFaces/2;

  CFreal faceCenterXPrevious = -1000.0, faceCenterXMin, faceCenterXCoordinate;

  // defining the connectivity between the periodic faces on the two TRs by sorting
  // the x-coordinates of their face centers from the smallest to the largest on each topological region  
  for (CFuint iFaceInSequence = 0; iFaceInSequence < nbTrFaces; ++iFaceInSequence) {
     faceCenterXMin = 1000.0;	
     for (CFuint iFace = 0; iFace < nbTrFaces; ++iFace) {
        // build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity *const face = geoBuilder->buildGE();
        const vector<Node*>& nodes = *face->getNodes();
	faceCenterXCoordinate = 0.5*((*(nodes[0]))[XX] + (*(nodes[1]))[XX]);
	cf_assert(faceCenterXCoordinate > -1000.0);
	cf_assert(faceCenterXCoordinate < 1000.0);
        if ((faceCenterXCoordinate < faceCenterXMin) && (faceCenterXCoordinate > faceCenterXPrevious)) {
       	  faceCenterXMin = faceCenterXCoordinate;
          _boundaryStatesInSequence[iFaceInSequence] = face->getState(0);
          const CFuint faceGlobalID = face->getID();
          _localTRSFaceIDInSequence.insert(faceGlobalID,iFaceInSequence);
     	}
     	// release the GeometricEntity
     	geoBuilder->releaseGE();
     }
     faceCenterXPrevious = faceCenterXMin;
  }

  faceCenterXPrevious = -1000.0;
  for (CFuint iFaceInSequence = nbTrFaces; iFaceInSequence < nbTrsFaces; ++iFaceInSequence) {
     faceCenterXMin = 1000.0;
     for (CFuint iFace = nbTrFaces; iFace < nbTrsFaces; ++iFace) {
        // build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity *const face = geoBuilder->buildGE();
        const vector<Node*>& nodes = *face->getNodes();
	faceCenterXCoordinate = 0.5*((*(nodes[0]))[XX] + (*(nodes[1]))[XX]);
        cf_assert(faceCenterXCoordinate > -1000.0);
        cf_assert(faceCenterXCoordinate < 1000.0);
        if ((faceCenterXCoordinate < faceCenterXMin) && (faceCenterXCoordinate > faceCenterXPrevious)) {
          faceCenterXMin = faceCenterXCoordinate;
          _boundaryStatesInSequence[iFaceInSequence] = face->getState(0);
          const CFuint faceGlobalID = face->getID();
          _localTRSFaceIDInSequence.insert(faceGlobalID,iFaceInSequence);
        }
        // release the GeometricEntity
        geoBuilder->releaseGE();
     }
     faceCenterXPrevious = faceCenterXMin;
  }

  _localTRSFaceIDInSequence.sortKeys();

  for (CFuint iPeriodicFace = 0; iPeriodicFace < nbTrFaces; ++iPeriodicFace) {
      _periodicFaceID[iPeriodicFace] = nbTrFaces + iPeriodicFace;
      const CFuint otherTRPeriodicFace = nbTrFaces + iPeriodicFace;
      _periodicFaceID[otherTRPeriodicFace] = iPeriodicFace;
  }
}
	
//////////////////////////////////////////////////////////////////////////////

void ShiftedPeriodicX2D::setGhostState(GeometricEntity *const face)
{
   State *const ghostState = face->getState(1);

   const CFuint faceGlobalID = face->getID();
   const CFuint periodicFaceID = _periodicFaceID[_localTRSFaceIDInSequence.find(faceGlobalID)];
   
   // watch out NOT to use the operator=, because in that case 
   // the overloaded version operator=(State) would be used =>
   // also the coordinates (Node) would be set equal!!! 

   ghostState->copyData(*_boundaryStatesInSequence[periodicFaceID]);
} 

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
