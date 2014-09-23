#include "FiniteVolume/FiniteVolume.hh"
#include "PeriodicY2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicY2D, CellCenterFVMData, FiniteVolumeModule> periodicY2DFVMCCProvider("PeriodicY2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void PeriodicY2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
}

//////////////////////////////////////////////////////////////////////////////

PeriodicY2D::PeriodicY2D(const std::string& name) :
  FVMCC_BC(name),
  _periodicFaceID(),
  _boundaryStates(),
  _globalToLocalTRSFaceID()
{
  addConfigOptionsTo(this);

  _threshold = 1e-5;
  setParameter("Threshold",&_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicY2D::~PeriodicY2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicY2D::setup()
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

  CFreal nodeBiggerYCoordinate  = 0.;
  CFreal nodeSmallerYCoordinate = 0.;
  CFreal nodePeriodicBiggerYCoordinate = 0.;
  CFreal nodePeriodicSmallerYCoordinate = 0.;

  // defining the connectivity between the periodic faces on the two TRs by checking
  // the equality of y-coordinates of their nodes

  for (CFuint iFace = 0; iFace < nbTrFaces; ++iFace) {
     // build the GeometricEntity
     geoData.idx = iFace;
     GeometricEntity *const face = geoBuilder->buildGE();
     const vector<Node*>& nodes = *face->getNodes();
     cf_assert(nodes.size() == 2);
     if((*(nodes[0]))[YY] > (*(nodes[1]))[YY]) {
       nodeBiggerYCoordinate = (*(nodes[0]))[YY];
       nodeSmallerYCoordinate = (*(nodes[1]))[YY];
     }
     else if((*(nodes[1]))[YY] > (*(nodes[0]))[YY]) {
	    nodeBiggerYCoordinate = (*(nodes[1]))[YY];
	    nodeSmallerYCoordinate = (*(nodes[0]))[YY];
          }
     // release the GeometricEntity
     geoBuilder->releaseGE();

     for (CFuint iPeriodicFace = nbTrFaces; iPeriodicFace < nbTrsFaces; ++iPeriodicFace) {
     	// build the GeometricEntity
     	geoData.idx = iPeriodicFace;
     	GeometricEntity *const facePeriodic = geoBuilder->buildGE();
     	const vector<Node*>& facePeriodicNodes = *facePeriodic->getNodes();
	cf_assert(facePeriodicNodes.size() == 2);
     	if((*(facePeriodicNodes[0]))[YY] > (*(facePeriodicNodes[1]))[YY]) {
          nodePeriodicBiggerYCoordinate = (*(facePeriodicNodes[0]))[YY];
          nodePeriodicSmallerYCoordinate = (*(facePeriodicNodes[1]))[YY];
        }
        else if((*(facePeriodicNodes[1]))[YY] > (*(facePeriodicNodes[0]))[YY]) {
	       nodePeriodicBiggerYCoordinate = (*(facePeriodicNodes[1]))[YY];
	       nodePeriodicSmallerYCoordinate = (*(facePeriodicNodes[0]))[YY];
        }
        // release the GeometricEntity
        geoBuilder->releaseGE();
        if(MathTools::MathChecks::isEqualWithError(nodeBiggerYCoordinate,nodePeriodicBiggerYCoordinate, _threshold) &&
	   MathTools::MathChecks::isEqualWithError(nodeSmallerYCoordinate,nodePeriodicSmallerYCoordinate, _threshold)) {
          _periodicFaceID[iFace] = iPeriodicFace;
	  _periodicFaceID[iPeriodicFace] = iFace;
          break;
        }
     }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicY2D::setGhostState(GeometricEntity *const face)
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
