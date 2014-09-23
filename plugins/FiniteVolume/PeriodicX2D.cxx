#include "MathTools/MathChecks.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/PeriodicX2D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicX2D, CellCenterFVMData, FiniteVolumeModule> periodicX2DFVMCCProvider("PeriodicX2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void PeriodicX2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
}

//////////////////////////////////////////////////////////////////////////////

PeriodicX2D::PeriodicX2D(const std::string& name) :
  FVMCC_BC(name),
  _periodicFaceID(),
  _boundaryStates(),
  _globalToLocalTRSFaceID()
{
  addConfigOptionsTo(this);

  _threshold = 1e-4;
  setParameter("Threshold",&_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicX2D::~PeriodicX2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicX2D::setup()
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
  _boundaryStates.resize(nbTrsFaces);

  // mapping the global face IDs to local IDs in the TRS and
  // storing the corresponding boundary states in a temporary vector

  CFreal nodeBiggerXCoordinate = 0.;
  CFreal nodeSmallerXCoordinate = 0.;
  CFreal nodePeriodicBiggerXCoordinate = 0.;
  CFreal nodePeriodicSmallerXCoordinate = 0.;



  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
     CFLogDebugMed( "iFace = " << iFace << "\n");
     // build the GeometricEntity
     geoData.idx = iFace;

     GeometricEntity *const face = geoBuilder->buildGE();
     _boundaryStates[iFace] = face->getState(0);
     const CFuint faceGlobalID = face->getID();
     _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);
     
     // release the GeometricEntity
     geoBuilder->releaseGE();
  }

  _globalToLocalTRSFaceID.sortKeys();
  const CFuint nbTrFaces = nbTrsFaces/2;

  
  // defining the connectivity between the periodic faces on the two TRs by checking
  // the equality of x-coordinates of their nodes

  for (CFuint iFace = 0; iFace < nbTrFaces; ++iFace) {
    // build the GeometricEntity
    geoData.idx = iFace;
    GeometricEntity *const face = geoBuilder->buildGE();
    const vector<Node*>& nodes = *face->getNodes();
    cf_assert(nodes.size() == 2);
    if((*(nodes[0]))[XX] > (*(nodes[1]))[XX]) {
      nodeBiggerXCoordinate = (*(nodes[0]))[XX];
      nodeSmallerXCoordinate = (*(nodes[1]))[XX];
    }
    else if((*(nodes[1]))[XX] > (*(nodes[0]))[XX]) {
      nodeBiggerXCoordinate = (*(nodes[1]))[XX];
      nodeSmallerXCoordinate = (*(nodes[0]))[XX];
    }
    // release the GeometricEntity
    geoBuilder->releaseGE();
    
    for (CFuint iPeriodicFace = nbTrFaces; iPeriodicFace < nbTrsFaces; ++iPeriodicFace) {
      // build the GeometricEntity
      geoData.idx = iPeriodicFace;
      GeometricEntity *const facePeriodic = geoBuilder->buildGE();
      const vector<Node*>& facePeriodicNodes = *facePeriodic->getNodes();
      //      const CFuint nbNodesInFacePeriodic = facePeriodicNodes.size();
      cf_assert(facePeriodicNodes.size() == 2);
      if((*(facePeriodicNodes[0]))[XX] > (*(facePeriodicNodes[1]))[XX]) {
	nodePeriodicBiggerXCoordinate = (*(facePeriodicNodes[0]))[XX];
	nodePeriodicSmallerXCoordinate = (*(facePeriodicNodes[1]))[XX];
      }
      else if((*(facePeriodicNodes[1]))[XX] > (*(facePeriodicNodes[0]))[XX]) {
	 nodePeriodicBiggerXCoordinate = (*(facePeriodicNodes[1]))[XX];
	 nodePeriodicSmallerXCoordinate = (*(facePeriodicNodes[0]))[XX];
       }
       // release the GeometricEntity
       geoBuilder->releaseGE();
       if(MathTools::MathChecks::isEqualWithError(nodeBiggerXCoordinate,nodePeriodicBiggerXCoordinate, _threshold) &&
	  MathTools::MathChecks::isEqualWithError(nodeSmallerXCoordinate,nodePeriodicSmallerXCoordinate, _threshold)) {
	 _periodicFaceID[iFace] = iPeriodicFace;
	 _periodicFaceID[iPeriodicFace] = iFace;
	 break;
       }
    }
  }
  
  // for (CFuint iFace = 0; iFace < nbTrFaces; ++iFace) {
//     // build the GeometricEntity
//     geoData.idx = iFace;
//     GeometricEntity *const face = geoBuilder->buildGE();
//     const vector<Node*>& nodes = *face->getNodes();    
//     cout <<  max( (*nodes[0])[XX], (*nodes[1])[XX]) << " ";
    
//     // release the GeometricEntity
//     geoBuilder->releaseGE();
    
//     geoData.idx = _periodicFaceID[iFace];
//     GeometricEntity *const oface = geoBuilder->buildGE();
//     const vector<Node*>& onodes = *oface->getNodes(); 
//     cout << max( (*onodes[0])[XX], (*onodes[1])[XX]) << ", " << iFace << " <=> " << _periodicFaceID[iFace] << endl << endl;
    
//     // release the GeometricEntity
//     geoBuilder->releaseGE();
//   }
  
//   abort();
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicX2D::setGhostState(GeometricEntity *const face)
{
   State *const ghostState = face->getState(1);

   const CFuint faceGlobalID = face->getID();
   const CFuint periodicFaceID = _periodicFaceID[_globalToLocalTRSFaceID.find(faceGlobalID)];

   // watch out NOT to use the operator=, because in that case
   // the overloaded version operator=(State) would be used =>
   // also the coordinates (Node) would be set equal!!!

   ghostState->copyData(*_boundaryStates[periodicFaceID]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
