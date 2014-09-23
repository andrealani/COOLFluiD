#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "UnsteadyMHD3DProjectionPhotosphere.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "Framework/DataHandle.hh"
#include "Framework/SubSystemStatus.hh"
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

MethodCommandProvider<UnsteadyMHD3DProjectionPhotosphere, CellCenterFVMData, FiniteVolumeMHDModule> 
unsteadyMHD3DProjectionPhotosphereFVMCCProvider("UnsteadyMHD3DProjectionPhotosphereFVMCC");

//////////////////////////////////////////////////////////////////////
   
void UnsteadyMHD3DProjectionPhotosphere::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("rhoFixed", "Non-dimensional density value that is to be fixed at the boundary.");
  options.addConfigOption< CFreal >("restartTime","Begin time of the overall unsteady simulation in case of restart (in seconds).");
  options.addConfigOption< CFreal >("initTime","Initial time of the overall unsteady simulation (in seconds).");
  options.addConfigOption< CFreal >("endTime","End time of the overall unsteady simulation (in seconds).");
}

//////////////////////////////////////////////////////////////////////

UnsteadyMHD3DProjectionPhotosphere::UnsteadyMHD3DProjectionPhotosphere(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  socket_BPFSS("BPFSS"),
  socket_BPFSSGhostBegin("BPFSSGhostBegin"),
  socket_BPFSSGhostEnd("BPFSSGhostEnd"),
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

  _restartTime = 0.0;
  setParameter("restartTime",&_restartTime);

  // this is always 0 as long as only two magnetogram data files are used per simulation
  _initTime = 0.0;
  setParameter("initTime",&_initTime);

  // 6 h = 21600 s, which is the frequency of the magnetogram data used in the calculation of the PFSS field
  // this is always 21600 as long as only two magnetogram data files are used per simulation
  _endTime = 21600.0;
  setParameter("endTime",&_endTime);

}

//////////////////////////////////////////////////////////////////////

UnsteadyMHD3DProjectionPhotosphere::~UnsteadyMHD3DProjectionPhotosphere()
{
}

//////////////////////////////////////////////////////////////////////

void UnsteadyMHD3DProjectionPhotosphere::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////

void UnsteadyMHD3DProjectionPhotosphere::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull());

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

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();

  CFout << "PFSS magnetic field is computed for " << trs->getName() << " on " << nbTrsFaces << " faces\n";

  DataHandle<std::vector<CFreal> > BPFSS = socket_BPFSS.getDataHandle();

  DataHandle<State*> gstates = socket_gstates.getDataHandle();
  const CFuint nbGhostStates = gstates.size();

  DataHandle<std::vector<CFreal> > BPFSSGhostBegin  = socket_BPFSSGhostBegin.getDataHandle();
  BPFSSGhostBegin.resize(nbGhostStates);
 
  DataHandle<std::vector<CFreal> > BPFSSGhostEnd  = socket_BPFSSGhostEnd.getDataHandle();
  BPFSSGhostEnd.resize(nbGhostStates);

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
                  geoBuilder = getMethodData().getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates,
                   socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  geoData.trs = trs;

  bool interpolationFlag;
  CFreal quadPointXCoord, quadPointYCoord, quadPointZCoord;
  // DataHandle does not accept RealVector so temporary vectors are created
  RealVector quadPointCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
          quadPointCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
	  ghostStateCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
          ghostStateCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
          BPFSSCartesian(PhysicalModelStack::getActive()->getDim()),
	  BPFSSGhostStateCartesian(PhysicalModelStack::getActive()->getDim());
  RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
          sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());

  vector<CFreal> BPFSSCartesianCoords(PhysicalModelStack::getActive()->getDim()),
	  BPFSSGhostStateCartesianCoords(PhysicalModelStack::getActive()->getDim());

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
          interpolationFlag = false;
	  _varSet->getModel()->setInterpolationFlag(interpolationFlag);
          _varSet->computePFSSMagneticField(quadPointCoordsSpherical,BPFSSCartesian,sphCarTransMat);

          BPFSSCartesianCoords[0] = BPFSSCartesian[0];
          BPFSSCartesianCoords[1] = BPFSSCartesian[1];
          BPFSSCartesianCoords[2] = BPFSSCartesian[2];

          BPFSS[faceID] = BPFSSCartesianCoords;

          // compute the PFSS magnetic field at the ghost state to be used for the interpolation

          State *const ghostState = face->getState(1);
	  const CFuint ghostStateID = ghostState->getLocalID();

	  ghostStateCoordsCartesian = ghostState->getCoordinates();

	  _varSet->setTransformationMatrices(ghostStateCoordsCartesian,ghostStateCoordsSpherical,carSphTransMat,sphCarTransMat);                    
          _varSet->computePFSSMagneticField(ghostStateCoordsSpherical,BPFSSGhostStateCartesian,sphCarTransMat);

          BPFSSGhostStateCartesianCoords[0] = BPFSSGhostStateCartesian[0];
	  BPFSSGhostStateCartesianCoords[1] = BPFSSGhostStateCartesian[1];
	  BPFSSGhostStateCartesianCoords[2] = BPFSSGhostStateCartesian[2]; 

	  BPFSSGhostBegin[ghostStateID] = BPFSSGhostStateCartesianCoords;

          interpolationFlag = true;
          _varSet->getModel()->setInterpolationFlag(interpolationFlag);
	  _varSet->computePFSSMagneticField(ghostStateCoordsSpherical,BPFSSGhostStateCartesian,sphCarTransMat);

	  BPFSSGhostStateCartesianCoords[0] = BPFSSGhostStateCartesian[0];
          BPFSSGhostStateCartesianCoords[1] = BPFSSGhostStateCartesian[1];
          BPFSSGhostStateCartesianCoords[2] = BPFSSGhostStateCartesian[2];

	  BPFSSGhostEnd[ghostStateID] = BPFSSGhostStateCartesianCoords;

          // release the GeometricEntity
          geoBuilder->releaseGE();
  }

}

//////////////////////////////////////////////////////////////////////

void UnsteadyMHD3DProjectionPhotosphere::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  DataHandle<std::vector<CFreal> > BPFSSGhostBegin  = socket_BPFSSGhostBegin.getDataHandle();
  DataHandle<std::vector<CFreal> > BPFSSGhostEnd  = socket_BPFSSGhostEnd.getDataHandle();

  const CFuint ghostStateID = ghostState->getLocalID();

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

  _dataGhostState[MHDProjectionTerm::RHO] = 2.*_rhoFixed - _dataInnerState[MHDProjectionTerm::RHO];
  _dataGhostState[MHDProjectionTerm::VX] = VCartesianGhostState[0];
  _dataGhostState[MHDProjectionTerm::VY] = VCartesianGhostState[1];
  _dataGhostState[MHDProjectionTerm::VZ] = VCartesianGhostState[2];
  _dataGhostState[MHDProjectionTerm::V] = VelCartesianGhostState;

  // interpolation to find BPFSS components at "time"

  const CFreal currentTime = SubSystemStatusStack::getActive()->getCurrentTime();
  const CFreal time = _restartTime+currentTime;

  cf_assert(time <= _endTime);
  cf_assert(time >= _initTime);
  cf_assert((BPFSSGhostBegin[ghostStateID])[0] != 0.0);
  cf_assert((BPFSSGhostBegin[ghostStateID])[1] != 0.0);
  cf_assert((BPFSSGhostBegin[ghostStateID])[2] != 0.0);
  cf_assert((BPFSSGhostEnd[ghostStateID])[0] != 0.0);
  cf_assert((BPFSSGhostEnd[ghostStateID])[1] != 0.0);
  cf_assert((BPFSSGhostEnd[ghostStateID])[2] != 0.0);
  cf_assert((BPFSSGhostBegin[ghostStateID])[0] != (BPFSSGhostEnd[ghostStateID])[0]);
  cf_assert((BPFSSGhostBegin[ghostStateID])[1] != (BPFSSGhostEnd[ghostStateID])[1]);
  cf_assert((BPFSSGhostBegin[ghostStateID])[2] != (BPFSSGhostEnd[ghostStateID])[2]);

  //cout << "TIME: " << time << " END_TIME: " << _endTime << " INIT_TIME: " << _initTime << endl;
  //cout << "BEGINX: " << (BPFSSGhostBegin[ghostStateID])[0] << " ENDX: " << (BPFSSGhostEnd[ghostStateID])[0] << endl;
  //cout << "BEGINY: " << (BPFSSGhostBegin[ghostStateID])[1] << " ENDY: " << (BPFSSGhostEnd[ghostStateID])[1] << endl;
  //cout << "BEGINZ: " << (BPFSSGhostBegin[ghostStateID])[2] << " ENDZ: " << (BPFSSGhostEnd[ghostStateID])[2] << endl << endl;
 
  const CFreal BxPFSSGhostTime = (BPFSSGhostBegin[ghostStateID])[0] + (((BPFSSGhostEnd[ghostStateID])[0]-(BPFSSGhostBegin[ghostStateID])[0])/(_endTime-_initTime))*(time-_initTime);  
  const CFreal ByPFSSGhostTime = (BPFSSGhostBegin[ghostStateID])[1] + (((BPFSSGhostEnd[ghostStateID])[1]-(BPFSSGhostBegin[ghostStateID])[1])/(_endTime-_initTime))*(time-_initTime);
  const CFreal BzPFSSGhostTime = (BPFSSGhostBegin[ghostStateID])[2] + (((BPFSSGhostEnd[ghostStateID])[2]-(BPFSSGhostBegin[ghostStateID])[2])/(_endTime-_initTime))*(time-_initTime);

  // assign the difference between the BPFSS components at "time" and "_initTime" to B1 at the ghost cell
  // !!! This is the driving mechanism of the unsteady simulations

  const CFreal Bx1GhostTime = BxPFSSGhostTime - (BPFSSGhostBegin[ghostStateID])[0];
  const CFreal By1GhostTime = ByPFSSGhostTime - (BPFSSGhostBegin[ghostStateID])[1];
  const CFreal Bz1GhostTime = BzPFSSGhostTime - (BPFSSGhostBegin[ghostStateID])[2];

  cf_assert(Bx1GhostTime == 0.0);
  cf_assert(By1GhostTime == 0.0);
  cf_assert(Bz1GhostTime == 0.0);

  _dataGhostState[MHDProjectionTerm::BX] = Bx1GhostTime;
  _dataGhostState[MHDProjectionTerm::BY] = By1GhostTime;
  _dataGhostState[MHDProjectionTerm::BZ] = Bz1GhostTime;
  const CFreal BGhostState = sqrt(_dataGhostState[MHDProjectionTerm::BX]*_dataGhostState[MHDProjectionTerm::BX] +
				_dataGhostState[MHDProjectionTerm::BY]*_dataGhostState[MHDProjectionTerm::BY] +
				_dataGhostState[MHDProjectionTerm::BZ]*_dataGhostState[MHDProjectionTerm::BZ]);
  _dataGhostState[MHDProjectionTerm::B] = BGhostState;

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

  const CFreal pGhostState = (nRef*k*2.0*TRef)/(rhoRef*vRef*vRef);

  _dataGhostState[MHDProjectionTerm::P] = pGhostState; 
  _dataGhostState[MHDProjectionTerm::A] = sqrt(_varSet->getModel()->getGamma()*_dataGhostState[MHDProjectionTerm::P]/_dataGhostState[MHDProjectionTerm::RHO]);
  _dataGhostState[MHDProjectionTerm::PHI] = _dataInnerState[MHDProjectionTerm::PHI];

  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
UnsteadyMHD3DProjectionPhotosphere::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
          FVMCC_BC::providesSockets();
  result.push_back(&socket_BPFSSGhostBegin);
  result.push_back(&socket_BPFSSGhostEnd);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > UnsteadyMHD3DProjectionPhotosphere::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();
  result.push_back(&socket_BPFSS);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
