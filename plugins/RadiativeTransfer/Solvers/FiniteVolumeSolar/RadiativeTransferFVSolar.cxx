#include <fstream>
#include <iostream>

#include "Common/PE.hh"
#include "Common/BadValueException.hh"
#include "Common/CFPrintContainer.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SocketBundleSetter.hh"

#include "FiniteVolume/CellCenterFVM.hh"

#include "RadiativeTransfer/RadiativeTransfer.hh"
#include "RadiativeTransfer/Solvers/FiniteVolumeSolar/RadiativeTransferFVSolar.hh"
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferFVSolar, DataProcessingData, RadiativeTransferModule>
radiativeTransferFVSolarProvider("RadiativeTransferFVSolar");

//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVSolar::RadiativeTransferFVSolar(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_faceCenters("faceCenters"),
  socket_faceAreas("faceAreas"),
  socket_divq("divq"),
  socket_qx("qx"),
  socket_qy("qy"),
  socket_qz("qz"),
  socket_qradFluxWall("qradFluxWall"),
  m_library(CFNULL),
  m_radiation(new RadiationPhysicsHandler("RadiationPhysicsHandler")),
  m_mapGeoToTrs(CFNULL),
  m_isWallFace(),
  m_cellID2WallFaceID(),
  m_tableWallFaceID2CellIDs(),
  m_mapWallTRSToOffsetFaceID(),
  m_cellBuilder(), 
  m_wallFaceBuilder(),
  m_normal()
{
  addConfigOptionsTo(this);
  
  /*m_maxNbNormalFaces = 10000;
    setParameter("maxNbNormalFaces",&m_maxNbNormalFaces); */
  
  // AL: to be removed once a cleaner solution is found 
  m_radNamespace = "Default";
  setParameter("RadNamespace", &m_radNamespace);
  
  m_wallTrsNames = vector<std::string>();
  setParameter("wallTrsNames",&m_wallTrsNames);
}
    
//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVSolar::~RadiativeTransferFVSolar()
{
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVSolar::defineConfigOptions(Config::OptionList& options)
{  
  // AL: to be removed once a cleaner solution is found 
  options.addConfigOption< string >
    ("RadNamespace","Namespace grouping all ranks involved in parallel communication");

  options.addConfigOption< vector<std::string> >
    ("wallTrsNames","Names of the TRSs from which the integration lines are constructed (solar surface)");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVSolar::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "RadiativeTransferFVSolar::configure() => START\n");
  
  DataProcessingCom::configure(args);
  
  cf_assert(m_radiation.isNotNull());
  configureNested ( m_radiation.getPtr(), args );
  
  CFLog(VERBOSE, "RadiativeTransferFVSolar::configure() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
RadiativeTransferFVSolar::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_divq);
  result.push_back(&socket_qx);
  result.push_back(&socket_qy);
  result.push_back(&socket_qz);
  result.push_back(&socket_qradFluxWall);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
RadiativeTransferFVSolar::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);  
  result.push_back(&socket_faceCenters);
  result.push_back(&socket_faceAreas);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVSolar::setup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "RadiativeTransferFVSolar::setup() => START\n");
  
  DataProcessingCom::setup();
  
  SocketBundle sockets;
  sockets.states      = socket_states;
  sockets.gstates     = socket_gstates;
  sockets.nstates     = socket_nstates;
  sockets.volumes     = socket_volumes;
  sockets.nodes       = socket_nodes;
  sockets.isOutward   = socket_isOutward;
  sockets.normals     = socket_normals;
  sockets.faceCenters = socket_faceCenters;
  sockets.faceAreas   = socket_faceAreas;
  
  m_radiation->setupDataSockets(sockets);
  m_radiation->setup();
  m_radiation->configureTRS();
  m_radiation->setupAxiFlag(false);
    
  // copy the content of the "initialNodalStates" (if available) into "nstates"
  const string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  const string dhName = nsp + "_initialNodalStates";
  if (MeshDataStack::getActive()->getDataStorage()->checkData(dhName)) {
    DataHandle<CFreal> initialNodalStates =
      MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(dhName); 
  
    /// storage of nstates (states in nodes)
    DataHandle<RealVector> nstates = socket_nstates.getDataHandle();  
    const CFuint nbNodalVars = nstates[0].size();
    cf_assert(nstates.size()*nbNodalVars == initialNodalStates.size());
    for (CFuint i = 0; i < nstates.size(); ++i) {
      const CFuint startn = i*nbNodalVars;
      for (CFuint j = 0; j < nbNodalVars; ++j) {
	nstates[i][j] = initialNodalStates[startn+j]; 
      }
    }  
    
    MeshDataStack::getActive()->getDataStorage()->deleteData<CFreal>(dhName);
  }
  
  /* if (m_radiation->hasRadiationPhysics()) {
    SafePtr<Radiator> radiator = 
      m_radiation->getCellDistPtr(0)->getRadiatorPtr();
    if (radiator->getName() != "NullRadiator") {
    }
   }
  */
  
  // cell builder
  m_cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_cellBuilder.setup();
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  cellData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  
  // wall face builder
  m_wallFaceBuilder.setup();
  m_wallFaceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_wallFaceBuilder.getDataGE().isBFace = true;
  
  m_mapGeoToTrs = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  
  CFLog(VERBOSE, "RadiativeTransferFVSolar::setup() => start\n");
  cf_assert(PE::GetPE().GetProcessorCount(nsp) == 1);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_normal.resize(dim, 0.); 
  
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();

  const CFuint nbCells = socket_states.getDataHandle().size(); 
  divQ.resize(nbCells);
  qx.resize(nbCells);
  qy.resize(nbCells);
  qz.resize(nbCells);
  
  // this is fundamental since += is used afterwards
  divQ   = 0.0;
  qx     = 0.0;
  qy     = 0.0;
  qz     = 0.0;
  
  // for which radiative heat flux has to be computed 
  const CFuint totalNbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  cf_assert(socket_normals.getDataHandle().size()/dim);
  m_isWallFace.resize(totalNbFaces, false);
  m_cellID2WallFaceID.resize(nbCells);
  
  CFuint nbFaces = 0; // total number of boundary faces belonging to TRS of type "Wall"  
  for(CFuint j=0; j< m_wallTrsNames.size(); ++j) {
    const string wallTRSName = m_wallTrsNames[j];
    SafePtr<TopologicalRegionSet> wallFaces = MeshDataStack::getActive()->getTrs(wallTRSName);
    CFLog(VERBOSE, "RadiativeTransferFVSolar::setup() => TRS["<< wallTRSName << "] is Wall\n");
    const CFuint nbLocalFaces = wallFaces->getLocalNbGeoEnts();
    for (CFuint f = 0; f < nbLocalFaces; ++f) {
      m_isWallFace[wallFaces->getLocalGeoID(f)] = true;
    }
    // store the offset for the corresponding wall TRS name
    m_mapWallTRSToOffsetFaceID[wallTRSName] = nbFaces;
    nbFaces += nbLocalFaces;
  }
  // preallocation of memory for qradFluxWall
  socket_qradFluxWall.getDataHandle().resize(nbFaces);
  socket_qradFluxWall.getDataHandle() = 0.0;
  
  Stopwatch<WallTime> stp;
  
  stp.start();
  storeIntegralPathIDs();
  CFLog(INFO, "RadiativeTransferFVSolar::setup() => getDirections() took " << stp.read() << "s\n");
  
  CFLog(VERBOSE, "RadiativeTransferFVSolar::setup() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVSolar::unsetup()
{
  CFAUTOTRACE;
  
  DataProcessingCom::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVSolar::storeIntegralPathIDs()
{
  CFLog(INFO, "RadiativeTransferFVSolar::storeIntegralPathIDs() => START\n");
  
  CFLog(VERBOSE, "RadiativeTransferFVSolar::storeIntegralPathIDs() => m_wallTrsNames.size() "
	<< m_wallTrsNames.size() << "\n");
  
  // this can only work serial!!
  if (m_wallTrsNames.size() > 0) {
    // locally built geo builders. At this stage there could be 
    // something missing that doesn't let use the ones owned by the method data
    // it is safer to use local ones
    GeometricEntityPool<FaceTrsGeoBuilder> faceBuilder;
    faceBuilder.setup();
    faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
    FaceTrsGeoBuilder::GeoData& faceData = faceBuilder.getDataGE();
    CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
    SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
    
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    // supposing that the mesh is fully cartesian/structured, we need mappings:
    // [cellID     -> wallfaceID] 
    // [wallfaceID -> cellIDs]
    CFuint nbWallFaces = 0;
    for (CFuint j = 0; j < m_wallTrsNames.size(); ++j) {
      SafePtr<TopologicalRegionSet> wallTRS = MeshDataStack::getActive()->getTrs(m_wallTrsNames[j]);
      nbWallFaces += wallTRS->getLocalNbGeoEnts(); 
    }
    
    // wall TRS ID (ID corresponding to the given wall face within the wall TRS list) -> faceID 
    vector<CFuint> wallTrsFaceID2FaceID(nbWallFaces);
    valarray<CFuint> nbCellsIDsPerWallTrsFaceID(nbWallFaces);
    nbCellsIDsPerWallTrsFaceID = 0;
    cf_assert(nbCellsIDsPerWallTrsFaceID.size() > 0);
    
    // first loop just to count and compute mapping [wallfaceID -> cellIDs]
    CFuint wallTrsFaceID = 0;
    for (CFuint j = 0; j < m_wallTrsNames.size(); ++j) {
      SafePtr<TopologicalRegionSet> wallTRS = MeshDataStack::getActive()->getTrs(m_wallTrsNames[j]);
      faceData.trs = wallTRS;
      
      CFuint cellID = 0;
      const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts(); 
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++wallTrsFaceID) {   
	faceData.isBFace = true;
	faceData.idx = iFace;
	GeometricEntity *const face = faceBuilder.buildGE();
	cellID = face->getState(0)->getLocalID();
	
	cf_assert(face->getState(1)->isGhost());
	CFuint faceID = face->getID();
	const CFuint wallFaceID = faceID;
	m_cellID2WallFaceID[cellID] = wallFaceID;
	nbCellsIDsPerWallTrsFaceID[wallTrsFaceID]++;
	wallTrsFaceID2FaceID[wallTrsFaceID] = wallFaceID;
	
	CFuint countFaces = 1;
	bool lastFace = false;
	while (!lastFace /*&& (countFaces < m_maxNbNormalFaces)*/) {
	  cellData.idx = cellID;
	  CFuint oppositeIFace = 0;
	  GeometricEntity *const cell = m_cellBuilder.buildGE();
	  const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos();
	  const CFuint nbCellFaces =  cellFaces.size();
	  
	  for (CFuint f = 0; f < nbCellFaces; ++f) {
	    if (cellFaces[f]->getID() == faceID) {
	      oppositeIFace = getOppositeIFace(f, PhysicalModelStack::getActive()->getDim(), 
					       cell->nbNodes());
	      countFaces++;
	      faceID = cellFaces[oppositeIFace]->getID();
	      
	      const CFuint leftID = cellFaces[oppositeIFace]->getState(0)->getLocalID();
	      if (!cellFaces[oppositeIFace]->getState(1)->isGhost()) {
		const CFuint rightID = cellFaces[oppositeIFace]->getState(1)->getLocalID();
		cellID = (leftID == cellID) ? rightID : leftID;
		nbCellsIDsPerWallTrsFaceID[wallTrsFaceID]++;
	      }
	      else {
		// if you reach a face with a ghost state, you have reached the opposite boundary: stop!
		lastFace = true;
		
		// max number of faces is set to the total number of faces
		/*m_maxNbNormalFaces = countFaces;*/
	      }
	      break;
	    }
	  }
	  m_cellBuilder.releaseGE();
	}
	faceBuilder.releaseGE();
	faceData.isBFace = false;
      }
    }
    
    if (nbCellsIDsPerWallTrsFaceID.sum() != nbCells) {
      CFLog(INFO, "nbCellsIDsPerWallTrsFaceID.sum() [" << nbCellsIDsPerWallTrsFaceID.sum() <<  "] != nbCells [" << nbCells << "]\n");
      cf_assert(nbCellsIDsPerWallTrsFaceID.sum() == nbCells);
    }
    
    m_tableWallFaceID2CellIDs.resize(nbCellsIDsPerWallTrsFaceID);
    cf_assert(m_tableWallFaceID2CellIDs.getSumCols() == nbCells);
    
    // ofstream fout("blFaces.dat");
    
    wallTrsFaceID = 0;
    for (CFuint j = 0; j < m_wallTrsNames.size(); ++j) {
      SafePtr<TopologicalRegionSet> wallTRS = 
	MeshDataStack::getActive()->getTrs(m_wallTrsNames[j]);
      faceData.trs = wallTRS;
      
      CFuint cellID = 0;
      const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts(); 
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++wallTrsFaceID) {   
	faceData.isBFace = true;
	faceData.idx = iFace;
	GeometricEntity *const face = faceBuilder.buildGE();
	cellID = face->getState(0)->getLocalID();
	
	cf_assert(face->getState(1)->isGhost());
	CFuint faceID = face->getID();
	const CFuint wallFaceID = faceID;
	CFuint iCell = 0;
	m_tableWallFaceID2CellIDs(wallTrsFaceID, iCell) = cellID;
	
	CFuint countFaces = 1;
	bool lastFace = false;
	while (!lastFace /*&& (countFaces < m_maxNbNormalFaces)*/) {
	  cellData.idx = cellID;
	  CFuint oppositeIFace = 0;
	  GeometricEntity *const cell = m_cellBuilder.buildGE();
	  const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos();
	  const CFuint nbCellFaces =  cellFaces.size();
	  
	  for (CFuint f = 0; f < nbCellFaces; ++f) {
	    if (cellFaces[f]->getID() == faceID) {
	      oppositeIFace = getOppositeIFace(f, PhysicalModelStack::getActive()->getDim(), 
					       cell->nbNodes());
	      countFaces++;
	      faceID = cellFaces[oppositeIFace]->getID();
	      
	      const CFuint leftID = cellFaces[oppositeIFace]->getState(0)->getLocalID();
	      if (!cellFaces[oppositeIFace]->getState(1)->isGhost()) {
		const CFuint rightID = cellFaces[oppositeIFace]->getState(1)->getLocalID();
		cellID = (leftID == cellID) ? rightID : leftID;
		m_tableWallFaceID2CellIDs(wallTrsFaceID, ++iCell) = cellID;
	      }
	      else {
		// if you reach a face with a ghost state, you have reached the opposite boundary: stop!
		lastFace = true;
		
		// max number of faces is set to the total number of faces
		/*m_maxNbNormalFaces = countFaces;*/
	      }
	      
	      // fout << cell->getState(0)->getCoordinates() << endl;
	      break;
	    }
	  }
	  m_cellBuilder.releaseGE();
	}
	// fout << endl;
	faceBuilder.releaseGE();
	faceData.isBFace = false;
      } 
    }
    // fout.close();
  }
  
  CFLog(INFO, "RadiativeTransferFVSolar::storeIntegralPathIDs() => END\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void RadiativeTransferFVSolar::execute()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "RadiativeTransferFVSolar::execute() => START\n");

  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  
  Stopwatch<WallTime> stp;
  stp.start();
  
  // this can only work serial!!
  if (m_wallTrsNames.size() > 0) {
    
    
    CFuint wallTrsFaceID = 0;
    // loop over photosphere TRS (typically is only one)
    for (CFuint j = 0; j < m_wallTrsNames.size(); ++j) {
      SafePtr<TopologicalRegionSet> wallTRS = MeshDataStack::getActive()->getTrs(m_wallTrsNames[j]);
      const CFuint nbTrsFaces = wallTRS->getLocalNbGeoEnts(); 
      CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
      
      // loop over photosphere TRS boundary faces
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++wallTrsFaceID) {
	const CFuint nbCellsInLine = m_tableWallFaceID2CellIDs.nbCols(wallTrsFaceID);
	for (CFuint iCell = 0; iCell < nbCellsInLine; ++iCell) {
	  const CFuint cellID = m_tableWallFaceID2CellIDs(wallTrsFaceID, iCell);
	  cellData.idx = cellID;
	  // "create" cell while marching normal to the solar photosphere
	  GeometricEntity *const cell = m_cellBuilder.buildGE();
	  const State& cellState = *cell->getState(0);

	  // if needed
	  const Node& coord = cellState.getCoordinates(); // position of cell center
	  const vector<Node*>& cellNodes = *cell->getNodes(); // vertices of cell 
	  
	  // here access density of whatevr you need from the state value
	  // ...
	  
	  // in case you need to use faces, you can access them this way
	  const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos();
	  //const CFuint nbCellFaces =  cellFaces.size();
	  //for (CFuint f = 0; f < nbCellFaces; ++f) {
	  //}
	  
	  // destroy cell
	  m_cellBuilder.releaseGE();
	}
      } 
    }
  }
  
  //reduceHeatFlux();
  
  
  CFLog(INFO, "RadiativeTransferFVSolar::execute() => took " << stp.read() << "s \n");
  CFLog(INFO, "RadiativeTransferFVSolar::execute() => END\n");
}
    
//////////////////////////////////////////////////////////////////////////////

CFuint RadiativeTransferFVSolar::getOppositeIFace(CFuint iFace, CFuint dim,
						  CFuint nbCellNodes) const
{
  if (dim == DIM_1D) {
   return (iFace == 0) ? 1 : 0;
  }
  else if (dim == DIM_2D && nbCellNodes == 4) {
    switch(iFace) {
    case 0:
      return 2;
      break;
    case 1:
      return 3;
      break;
    case 2:
      return 0;
      break;
    case 3:
      return 1;
      break;
    }
  }
  else if (dim == DIM_3D && nbCellNodes == 8) {
    switch(iFace) {
    case 0:
      return 1;
      break;
    case 1:
      return 0;
      break;
    case 2:
      return 4;
      break;
    case 3:
      return 5;
      break;
    case 4:
      return 2;
      break;
    case 5:
      return 3;
      break;
    }
  }
  else if (dim == DIM_3D && nbCellNodes == 6) {
    switch(iFace) {
    case 0:
      return 1;
      break;
    case 1:
      return 0;
      break;
    case 2:
      CFLog(INFO, "Prism face 2 doesn't have a single opposite\n");
      cf_assert(nbCellNodes == 6);
      break;
    case 3:
      CFLog(INFO, "Prism face 3 doesn't have a single opposite\n");
      cf_assert(nbCellNodes == 6);
      break;
    case 4:
      CFLog(INFO, "Prism face 4 doesn't have a single opposite\n");
      cf_assert(nbCellNodes == 6);
      break;
    }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

/*void RadiativeTransferFVSolar::getAdvanceOrder(const CFuint d, 
					     CFint *const advanceOrder)
{
  CFLog(VERBOSE, "RadiativeTransferFVSolar::getAdvanceOrder() => start\n");
  
  // The order of advance calculation begins here
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  const CFuint nbCells = CellID.size();
  cf_assert(nbCells > 0);
  const CFuint DIM = PhysicalModelStack::getActive()->getDim();
  cf_assert(DIM == DIM_3D);
  
  CFLog(INFO, "RadiativeTransferFVSolar::getAdvanceOrder() => Direction number [" << d <<"]\n");
  
  // precompute the dot products for all faces and directions (a part from the sign)
  computeDotProdInFace(d, m_dotProdInFace);
  
  CFuint mLast = 0;
  CFuint m = 0;

  CFuint stage = 1;
  // Initializing the sdone and cdone for each direction
  m_sdone.assign(nbCells, false);
  
  m_cdoneIdx.clear();
  cf_assert(m_cdoneIdx.size() == 0);
  cf_assert(m_cdoneIdx.capacity() == nbCells);
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces = MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
    
  while (m < nbCells) { //The loop over the cells begins
    mLast = m;	  //Checking to see if it counts all the cells
    for (CFuint iCell = 0; iCell < nbCells; iCell++) {
      CFLog(DEBUG_MAX, "RadiativeTransferFVSolar::getAdvanceOrder() => iCell = " << iCell <<"\n");
      if (m_sdone[iCell] == false) {
	const CFuint nbFaces = cellFaces->nbCols(iCell);
	for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	  const CFuint faceID = (*cellFaces)(iCell, iFace);
	  const bool isBFace  = m_mapGeoToTrs->isBGeo(faceID);
	  const bool neighborIsSdone = (!isBFace) ? m_sdone[getNeighborCellID(faceID, iCell)] : true;
	  
	  // When dot product is < 0 and neighbor is undone, skip the cell and continue the loop
	  const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	  const CFreal dotMult = m_dotProdInFace[faceID]*factor;
	  if (dotMult < 0. && neighborIsSdone == false) {goto cell_loop;}
	}// end loop over the FACES
	
	CFLog(DEBUG_MAX, "advanceOrder[" << d << "][" << m <<"] = " << iCell << "\n");
	
	advanceOrder[m] = iCell;
	CellID[iCell] = stage;
	m_cdoneIdx.push_back(iCell);
	m += 1;
      }// end if(Cell is not done)
      
    cell_loop:
      const bool dummy = true;
    }// end of the loop over the CELLS
    
    const string msg = "advanceOrder[" + StringOps::to_str(d) + "] = ";
    CFLog(DEBUG_MAX, msg << "\n");
    for (CFuint a = 0; a < nbCells; ++a) {
      CFLog(DEBUG_MAX, advanceOrder[a] << " ");
    }
    CFLog(DEBUG_MAX, "\n");

    if (m == mLast) {diagnoseProblem(d, m, mLast);}
    
    advanceOrder[m - 1] *= -1;
    
    CFLog(DEBUG_MAX, "At stage[" << stage << "] => m_cdoneIdx.size() = " << m_cdoneIdx.size() << "\n");
    for (CFuint id = 0; id < m_cdoneIdx.size(); ++id) {
      m_sdone[m_cdoneIdx[id]] = true;
    }
    m_cdoneIdx.clear();
    
    CFLog(VERBOSE, "RadiativeTransferFVSolar::getAdvanceOrder() => m  "<< m << " \n");
    CFLog(VERBOSE, "RadiativeTransferFVSolar::getAdvanceOrder() => End of the "<< stage << " stage\n");
    
    ++stage;
  }// end of the loop over the STAGES
  
  //Printing advanceOrder for debug purpuses
  const string msg = "RadiativeTransferFVSolar::getAdvanceOrder() => advanceOrder[" + StringOps::to_str(d) + "] = ";
  
  CFLog(DEBUG_MAX, msg << "\n");
  for (CFuint a = 0; a < nbCells; ++a) {
    CFLog(DEBUG_MAX, advanceOrder[a] << " ");
  }
  CFLog(DEBUG_MAX, "\n");
    
  CFLog(VERBOSE, "RadiativeTransferFVSolar::getAdvanceOrder() => end\n");
}*/
      
//////////////////////////////////////////////////////////////////////////////

/*void RadiativeTransferFVSolar::reduceHeatFlux()
{
  if (PE::GetPE().GetProcessorCount(m_radNamespace) > 1) {
    CFLog(VERBOSE, "RadiativeTransferFVSolar::reduceHeatFlux() => start\n");
    
    // AL: the following has to be moved to a postprocess() or something function
    //     to be called after execute() unless the ConcurrentCoupler can be modified to
    //     take care of it
    DataHandle< CFreal > qradFluxWall = socket_qradFluxWall.getDataHandle(); 
    const CFuint qsize = qradFluxWall.size();
    if (qsize > 0) { 
      PE::GetPE().setBarrier(m_radNamespace);
      
      vector<CFreal> localHeatFlux(qsize);
      for (CFuint i = 0; i < qsize; ++i) {
        localHeatFlux[i] = qradFluxWall[i];
      }
      
      MPIError::getInstance().check
	("MPI_Allreduce", "RadiativeTransferFVSolar::reduceHeatFlux()",
	 MPI_Allreduce(&localHeatFlux[0], &qradFluxWall[0], qsize, 
		       MPIStructDef::getMPIType(&localHeatFlux[0]), 
		       MPI_SUM, PE::GetPE().GetCommunicator(m_radNamespace)));
    }
    
    
    CFLog(VERBOSE, "RadiativeTransferFVSolar::reduceHeatFlux() => end\n");
  }
  }*/

//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

