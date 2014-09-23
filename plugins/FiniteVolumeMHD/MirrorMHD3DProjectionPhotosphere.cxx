#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MirrorMHD3DProjectionPhotosphere.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodCommandProvider<MirrorMHD3DProjectionPhotosphere, CellCenterFVMData, FiniteVolumeMHDModule> 
mirrorMHD3DProjectionPhotosphereFVMCCProvider("MirrorMHD3DProjectionPhotosphereFVMCC");

//////////////////////////////////////////////////////////////////////
   
void MirrorMHD3DProjectionPhotosphere::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("rhoFixed", "Non-dimensional density value that is to be fixed in the ghost cells.");
}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionPhotosphere::MirrorMHD3DProjectionPhotosphere(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  socket_BPFSS("BPFSS"),
  _dataInnerState(),
  _dataGhostState(),
  _cartesianSphericalTMInnerState(),
  _sphericalCartesianTMInnerState(),
  _cartesianSphericalTMGhostState(),
  _sphericalCartesianTMGhostState()
{
  addConfigOptionsTo(this);

  _rhoFixed = 1.0;
  setParameter("rhoFixed",&_rhoFixed);

}

//////////////////////////////////////////////////////////////////////

MirrorMHD3DProjectionPhotosphere::~MirrorMHD3DProjectionPhotosphere()
{
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPhotosphere::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPhotosphere::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull());

  const std::string potentialBType = _varSet->getPotentialBType();

  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _cartesianSphericalTMInnerState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
             Framework::PhysicalModelStack::getActive()->getDim());
  _sphericalCartesianTMInnerState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
                Framework::PhysicalModelStack::getActive()->getDim());
  _cartesianSphericalTMGhostState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
             Framework::PhysicalModelStack::getActive()->getDim());
  _sphericalCartesianTMGhostState.resize(Framework::PhysicalModelStack::getActive()->getDim(),
                Framework::PhysicalModelStack::getActive()->getDim());

  if (potentialBType == "PFSS") {

        SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
        const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();

        bool interpolationFlag = false;
        _varSet->getModel()->setInterpolationFlag(interpolationFlag);

        CFout << "PFSS magnetic field is computed for " << trs->getName() << " on " << nbTrsFaces << " faces\n";

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

//////////////////////////////////////////////////////////////////////

void MirrorMHD3DProjectionPhotosphere::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  RealVector innerStateCoordsSpherical(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector ghostStateCoordsSpherical(Framework::PhysicalModelStack::getActive()->getDim());

  const RealVector innerStateCoords = innerState->getCoordinates();
  const RealVector ghostStateCoords = ghostState->getCoordinates();

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  // set the transformation matrices between Cartesian and spherical coordinate systems
  _varSet->setTransformationMatrices(innerStateCoords,innerStateCoordsSpherical,_cartesianSphericalTMInnerState,_sphericalCartesianTMInnerState);
  _varSet->setTransformationMatrices(ghostStateCoords,ghostStateCoordsSpherical,_cartesianSphericalTMGhostState,_sphericalCartesianTMGhostState);

  RealVector VCartesianInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VSphericalInnerState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VCartesianGhostState(Framework::PhysicalModelStack::getActive()->getDim());
  RealVector VSphericalGhostState(Framework::PhysicalModelStack::getActive()->getDim());

  VCartesianInnerState[0] = _dataInnerState[MHDProjectionTerm::VX];
  VCartesianInnerState[1] = _dataInnerState[MHDProjectionTerm::VY];
  VCartesianInnerState[2] = _dataInnerState[MHDProjectionTerm::VZ];

  VSphericalInnerState = _cartesianSphericalTMInnerState*VCartesianInnerState;

  VSphericalGhostState[0] = 0.0;
  VSphericalGhostState[1] = VSphericalInnerState[1];
  VSphericalGhostState[2] = -VSphericalInnerState[2];

  VCartesianGhostState = _sphericalCartesianTMGhostState*VSphericalGhostState;

  const CFreal VelCartesianGhostState = sqrt(VCartesianGhostState[0]*VCartesianGhostState[0]
			      + VCartesianGhostState[1]*VCartesianGhostState[1]
			      + VCartesianGhostState[2]*VCartesianGhostState[2]);

  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  const CFreal nRef = _varSet->getNRef();
  const CFreal TRef = _varSet->getTRef();

  const CFreal rhoRef = nRef*(mp+me);
  const CFreal vRef = sqrt(2.0*k*TRef/mp);

  const CFreal pGhostState = (nRef*k*TRef)/(rhoRef*vRef*vRef); 

  _dataGhostState[MHDProjectionTerm::RHO] = 2.*_rhoFixed - _dataInnerState[MHDProjectionTerm::RHO];
  _dataGhostState[MHDProjectionTerm::VX] = VCartesianGhostState[0];
  _dataGhostState[MHDProjectionTerm::VY] = VCartesianGhostState[1];
  _dataGhostState[MHDProjectionTerm::VZ] = VCartesianGhostState[2];
  //!!! FOR THE MOMENT, THE VARIABLE MAGNETIC FIELD IS SIMPLY ASSIGNED TO BE ZERO SO THAT B=B0 AT THE BOUNDARY CELL FACE
  _dataGhostState[MHDProjectionTerm::BX] = -_dataInnerState[MHDProjectionTerm::BX]; 
  _dataGhostState[MHDProjectionTerm::BY] = -_dataInnerState[MHDProjectionTerm::BY];
  _dataGhostState[MHDProjectionTerm::BZ] = -_dataInnerState[MHDProjectionTerm::BZ];
  _dataGhostState[MHDProjectionTerm::V] = VelCartesianGhostState;
  _dataGhostState[MHDProjectionTerm::P] = pGhostState;
  _dataGhostState[MHDProjectionTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDProjectionTerm::P]/_dataGhostState[MHDProjectionTerm::RHO]);
  _dataGhostState[MHDProjectionTerm::B] = _dataInnerState[MHDProjectionTerm::B]; 
  _dataGhostState[MHDProjectionTerm::PHI] = _dataInnerState[MHDProjectionTerm::PHI];
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > MirrorMHD3DProjectionPhotosphere::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();
  result.push_back(&socket_BPFSS);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
