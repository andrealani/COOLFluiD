#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "SuperOutletMHD3DProjectionHeating.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperOutletMHD3DProjectionHeating, CellCenterFVMData, FiniteVolumeMHDModule> 
superOutletMHD3DProjectionHeatingFVMCCProvider("SuperOutletMHD3DProjectionHeatingFVMCC");
    
//////////////////////////////////////////////////////////////////////////////
   
void SuperOutletMHD3DProjectionHeating::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD3DProjectionHeating::SuperOutletMHD3DProjectionHeating(const std::string& name) : 
  FVMCC_BC(name),
  _varSet(CFNULL),
  socket_BPFSS("BPFSS"),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD3DProjectionHeating::~SuperOutletMHD3DProjectionHeating() 
{
}

//////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjectionHeating::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull());

  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  const std::string potentialBType = _varSet->getPotentialBType();

  if (potentialBType == "PFSS") {

  	SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  	const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();

  	CFout << "PFSS magnetic field is computed for " << trs->getName() << " on " << nbTrsFaces << " faces\n";

  	bool interpolationFlag = false;
  	_varSet->getModel()->setInterpolationFlag(interpolationFlag);

  	DataHandle<std::vector<CFreal> > BPFSS = socket_BPFSS.getDataHandle();

  	Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
         		geoBuilder = getMethodData().getFaceTrsGeoBuilder();

  	SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  	geoBuilderPtr->setDataSockets(socket_states, socket_gstates,
             	               socket_nodes);

  	FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  	geoData.isBFace = true;
  	geoData.trs = trs;

  	CFreal quadPointXCoord, quadPointYCoord, quadPointZCoord;
  	// DataHandle does not accept RealVector so temporary vectors are created
  	RealVector quadPointCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
        	quadPointCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
        	BPFSSCartesian(PhysicalModelStack::getActive()->getDim());
  	RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
        	sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());

  	vector<CFreal> BPFSSCartesianCoords(PhysicalModelStack::getActive()->getDim());

  	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
   		CFLogDebugMed( "iFace = " << iFace << "\n");

   		// build the GeometricEntity
   		geoData.idx = iFace;

   		GeometricEntity *const face = geoBuilder->buildGE();
   		const vector<Node*>& nodes = *face->getNodes();

   		const CFuint nbNodesInFace = nodes.size();

   		quadPointXCoord = 0.0;
   		quadPointYCoord = 0.0;
   		quadPointZCoord = 0.0;

   		for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
        		quadPointXCoord += (*(nodes[iNode]))[XX];
        		quadPointYCoord += (*(nodes[iNode]))[YY];
        		quadPointZCoord += (*(nodes[iNode]))[ZZ];
   		}	

   		quadPointXCoord /= nbNodesInFace;
   		quadPointYCoord /= nbNodesInFace;
   		quadPointZCoord /= nbNodesInFace;

   		quadPointCoordsCartesian[0] = quadPointXCoord;
   		quadPointCoordsCartesian[1] = quadPointYCoord;
   		quadPointCoordsCartesian[2] = quadPointZCoord;

   		const CFuint faceID = face->getID();

   		_varSet->setTransformationMatrices(quadPointCoordsCartesian,quadPointCoordsSpherical,carSphTransMat,sphCarTransMat);
   		_varSet->computePFSSMagneticField(quadPointCoordsSpherical,BPFSSCartesian,sphCarTransMat);

   		BPFSSCartesianCoords[0] = BPFSSCartesian[0];
   		BPFSSCartesianCoords[1] = BPFSSCartesian[1];
   		BPFSSCartesianCoords[2] = BPFSSCartesian[2];

   		BPFSS[faceID] = BPFSSCartesianCoords;

   		// release the GeometricEntity
   		geoBuilder->releaseGE();
  	}

  }
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjectionHeating::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD3DProjectionHeating::setGhostState(GeometricEntity *const face)
 {
   State *const innerState = face->getState(0);
   State *const ghostState = face->getState(1);
   // watch out NOT to use the operator=, because in that case 
   // the overloaded version operator=(State) would be used =>
   // also the coordinates (Node) would be set equal!!! 
   
   // set the physical data starting from the inner state
   _varSet->computePhysicalData(*innerState, _dataInnerState);

   _dataGhostState = _dataInnerState;

   // there are two possible boundary conditions for phi

   // 1) a reference value should be imposed for phi
   _dataGhostState[MHDProjectionTerm::PHI] = _refPhi;

   // 2)
   //_dataGhostState[MHDProjectionTerm::PHI] = -_dataInnerState[MHDProjectionTerm::PHI];

   _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > SuperOutletMHD3DProjectionHeating::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();
  result.push_back(&socket_BPFSS);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
