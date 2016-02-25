#include "MathTools/RealMatrix.hh"
#include "MathTools/InverterT.hh"

#include "Framework/MeshData.hh"
#include "Framework/ComputeNormals.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/VolumeCalculator.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/ComputeFaceNormalsFVMCC.hh"
#include "Framework/ComputeDummyStates.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/ElementDataArray.hh"
#include "Framework/SubSystemStatus.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/ComputeStencil.hh"

#include "Common/PE.hh"
#include "Common/CFMultiMap.hh"

#ifdef CF_HAVE_MPI
#include <mpi.h>
#include "Common/MPI/MPIStructDef.hh"
#endif

#include "FiniteVolumeCUDA/StencilCUDASetup.hh"
#include "FiniteVolumeCUDA/FiniteVolumeCUDA.hh"

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

MethodCommandProvider<StencilCUDASetup, CellCenterFVMData,FiniteVolumeCUDAModule> 
stencipCUDASetupProvider("StencilCUDASetup");

//////////////////////////////////////////////////////////////////////////////

void StencilCUDASetup::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("stencil","stencil type");
  options.addConfigOption< std::string >
    ("InitLimiterSocket","name of the input limiter socket");
}

//////////////////////////////////////////////////////////////////////////////

StencilCUDASetup::StencilCUDASetup(const std::string& name) :
  CellCenterFVMCom(name),
  _dynamicSockets(),
  socket_isOutward("isOutward"),
  socket_cellFlag("cellFlag"),
  socket_normals("normals"),
  socket_faceAreas("faceAreas"),  
  socket_gstates("gstates"),
  socket_volumes("volumes"),
  socket_limiter("limiter"),
  socket_nstates("nstates"),
  socket_nstatesProxy("nstatesProxy"),
  socket_trsID("trsID"),
  socket_rankPartitionFaces("rankPartitionFaces"),
  socket_stencil("stencil"),
  socket_activeDiffusion("activeDiffusion"),
  socket_states("states"),
  socket_nodes("nodes"),
  m_stored_args()
{
  addConfigOptionsTo(this);
  
  _stencilType = "FaceVertexPlusGhost";
  setParameter("stencil",&_stencilType);
  
  _limiterSocketName = "Null";
  setParameter("InitLimiterSocket",&_limiterSocketName);
}

//////////////////////////////////////////////////////////////////////////////

StencilCUDASetup::~StencilCUDASetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > StencilCUDASetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_cellFlag);
  result.push_back(&socket_gstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_limiter);
  result.push_back(&socket_nstates);
  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_trsID);
  result.push_back(&socket_rankPartitionFaces);
  result.push_back(&socket_stencil);
  result.push_back(&socket_activeDiffusion);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StencilCUDASetup::configure ( Config::ConfigArgs& args )
{
  CellCenterFVMCom::configure(args);

  //Loop over the TRS's and add the "TRSName" + "-boundaryNormals" datasocketsink to the _dynamicSockets
  
  std::string name = getMethodData().getNamespace();
  const string ssName = SubSystemStatusStack::getCurrentName();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance(ssName).getNamespace(name);
  Common::SafePtr<MeshData> meshData = MeshDataStack::getInstance().getEntryByNamespace(nsp);
  
  vector<std::string> trsList = meshData->getTRSNameList();
  const CFint nbTRSs = trsList.size();
  
  for (CFint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    const std::string trsName = trsList[iTRS];
    
    if (trsName != "PartitionFaces" &&
        trsName != "InnerCells" &&
        trsName != "InnerFaces")
      {
	const std::string socketName = trsName + "-boundaryNormals";
	const bool isEssential = false;
	_dynamicSockets.createSocketSink<const CFreal*>(socketName,isEssential);
      }
  }
  
  // here call the mesh reader to ask if _limiterSocketName is one of the extra fields
  //  getMethodData().getCollaborator<MeshCreator>()->getMethodData();
  if (_limiterSocketName != "Null") {
    _dynamicSockets.createSocketSink<CFreal>(_limiterSocketName, true);
  }
  
  m_stored_args = args;
  
  // set the number of overlap layers to two for parallel runs
  if (PE::GetPE().IsParallel()) {
    meshData->setNbOverlapLayers(2);
  }
}      

//////////////////////////////////////////////////////////////////////////////

void StencilCUDASetup::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "StencilCUDASetup::execute() => start\n");
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  isOutward.resize(nbFaces);
  isOutward = -1;
  
  DataHandle<bool> cellFlag = socket_cellFlag.getDataHandle();
  cellFlag.resize(states.size());
  cellFlag = false;
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  normals.resize(nbFaces*PhysicalModelStack::getActive()->getDim());
  
  DataHandle<CFreal> faceAreas = socket_faceAreas.getDataHandle();
  faceAreas.resize(nbFaces);
  
  computeNormalsData();

  DataHandle<State*> gstates = socket_gstates.getDataHandle();

  // compute the dummy (ghost) states
  ComputeDummyStates computeDummyStates;
  computeDummyStates.setDataSockets(socket_normals, socket_gstates,socket_states, socket_nodes);
  computeDummyStates.setup();
  computeDummyStates(getMethodData().getTRSsWithGhostsOnFace());
  
  if (PE::GetPE().IsParallel ()) {
    CFLog(VERBOSE, "StencilCUDASetup::execute() => before assignPartitionFaceGlobalGhostStateIDS()\n");
    assignPartitionFaceGlobalGhostStateIDS();
    CFLog(VERBOSE, "StencilCUDASetup::execute() => end assignPartitionFaceGlobalGhostStateIDS()\n");
  }
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  volumes.resize(states.size());
  computeCellVolumes();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // it must have size == 0 in this case
  DataHandle<CFreal> limiter = socket_limiter.getDataHandle();
  const CFuint sizeLimiter = states.size()*nbEqs;
  limiter.resize(sizeLimiter);
  
  CFLog(VERBOSE, "_limiterSocketName " << _limiterSocketName << "\n");
  
  if (_dynamicSockets.sinkSocketExists(_limiterSocketName)) {
    DataHandle<CFreal> initLimiter = _dynamicSockets.
      getSocketSink<CFreal>(_limiterSocketName)->getDataHandle();
    
    cf_assert(initLimiter.size() == limiter.size());
    for (CFuint i = 0; i < sizeLimiter; ++i) {
      limiter[i] = initLimiter[i];
    }
  }
  
  DataHandle<CFreal> activeDiffusion = socket_activeDiffusion.getDataHandle();
  activeDiffusion.resize(socket_states.getDataHandle().size());
  activeDiffusion = 1.;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  nstates.resize(nodes.size());
  for (CFuint i = 0; i < nstates.size(); ++i) {
    nstates[i].resize(nbEqs);
  }

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();
  nstatesProxy.resize(1);
  nstatesProxy[0] = new DofDataHandleIterator<RealVector,RealVector>(nstates);
  
  DataHandle<CFint> trsID = socket_trsID.getDataHandle();
  trsID.resize(nodes.size());
  trsID = -1;
  
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  stencil.resize(states.size());

  // ghost states are FOR THE MOMENT excluded from the stencil
  computeStencil();
  
  CFLog(VERBOSE, "StencilCUDASetup::execute() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void StencilCUDASetup::computeNormalsData()
{
  CFAUTOTRACE;

  SafePtr<vector<ElementTypeData> > elemTypes =
    MeshDataStack::getActive()->getElementTypeData();

  DataSocketSink<CFreal> sinkNormals(socket_normals);
  DataSocketSink<CFint> sinkIsOutward(socket_isOutward);

  SafePtr<DataSocketSink<CFreal> > sinkNormalsPtr = &sinkNormals;
  SafePtr<DataSocketSink<CFint> > sinkIsOutwardPtr = &sinkIsOutward;

  for (CFuint iType = 0; iType < elemTypes->size(); ++iType)
  {

    const CFuint geoOrder = (*elemTypes)[iType].getGeoOrder();

    const std::string elemName = (*elemTypes)[iType].getShape() + CFPolyOrder::Convert::to_str(geoOrder);

    SelfRegistPtr<ComputeNormals> computeFaceNormals =
      Environment::Factory<ComputeNormals>::getInstance().
      getProvider("Face" + elemName)->create();

    const CFuint firstElem = (*elemTypes)[iType].getStartIdx();
    const CFuint lastElem  = (*elemTypes)[iType].getEndIdx();

    SelfRegistPtr<ComputeFaceNormalsFVMCC> faceNormalsComputer =
      computeFaceNormals.d_castTo<ComputeFaceNormalsFVMCC>();

    faceNormalsComputer->setSockets(sinkNormalsPtr,
				    sinkIsOutwardPtr);

    (*faceNormalsComputer)(firstElem, lastElem);
  }
  
  // computation of face areas
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  DataHandle< CFreal> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> faceAreas = socket_faceAreas.getDataHandle();
  RealVector faceNormal(dim);
  
  for (CFuint iFace = 0; iFace < faceAreas.size(); ++iFace) {
    const CFuint startID = iFace*dim;
    for (CFuint i = 0; i < dim; ++i) {      
      faceNormal[i] = normals[startID + i];
    }
    faceAreas[iFace] = faceNormal.norm2();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StencilCUDASetup::computeCellVolumes()
{
  CFAUTOTRACE;

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
   getTrs("InnerCells");

  Common::SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> >
    geoBuilder = getMethodData().getGeoWithNodesBuilder();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;
  
  CFuint countNegativeVol = 0;
  const CFuint nbElems = cells->getLocalNbGeoEnts();
  for (CFuint iElem = 0; iElem < nbElems; ++iElem) {
    // build the GeometricEntity
    geoData.idx = iElem;
    GeometricEntity *const cell = geoBuilder->buildGE();
    volumes[iElem] = cell->computeVolume();
    
    if (volumes[iElem] < 0.0)
    {
      countNegativeVol++;
      cout.precision(14);
      CFout << "Cell [" << iElem << "] with [" << cell->nbNodes() << "] nodes has negative volume [" << volumes[iElem] << "]\n";
      volumes[iElem] = 0.0;
    }

    if ( MathChecks::isZero(volumes[iElem]) )
    {
      cout.precision(14);
      CFout << "Cell [" << iElem << "] with [" << cell->nbNodes() << "] nodes has zero volume [" << volumes[iElem] << "]\n";
      // print coordinates
      for (CFuint i = 0; i < cell->nbNodes(); ++i) {
       CFout << *cell->getNode(i) << ", ";
      } CFout << "\n";
      volumes[iElem] = 0.0;
    }

    //   cf_assert(volumes[iElem] > 0.);

    //release the GeometricEntity
    geoBuilder->releaseGE();
  }
  
  CFLogDebugMin ( "Negative volumes count [" << countNegativeVol << "]\n" );
  /// @todo Instead of asserting, add a user option
  ///       to accept continuing with negative volumes
  //  cf_assert(countNegativeVol == 0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
StencilCUDASetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result  = _dynamicSockets.getAllSinkSockets();

  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
      
void StencilCUDASetup::assignPartitionFaceGlobalGhostStateIDS()
{  
#ifdef CF_HAVE_MPI
  
  GeometricEntityPool<FaceTrsGeoBuilder> faceBuilder;
  faceBuilder.setup();
  faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = faceBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> pTRS = MeshDataStack::getActive()->getTrs("PartitionFaces");
  faceData.trs = pTRS;
  faceData.isBFace = true;
  
  // count the entries into the soon-to-be-built map in order to save memory 
  const CFuint nbFaces = faceData.trs->getLocalNbGeoEnts();
  CFuint faceNodesSize = 0;
  for (CFuint i = 0; i < nbFaces; ++i) {
    faceNodesSize += pTRS->getNbNodesInGeo(i);
  }
  
  DataHandle<CFuint> rankPartitionFaces = socket_rankPartitionFaces.getDataHandle();
  rankPartitionFaces.resize(nbFaces);
  
  // create a mapping between partition face nodes and corresponding local TRS face ID
  CFMultiMap<CFuint, CFuint> mapNodeIDToTrsFaceID(faceNodesSize);
  typedef CFMultiMap<CFuint, CFuint>::MapIterator MapIt;
  
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    faceData.idx = iFace;
    GeometricEntity *const face = faceBuilder.buildGE();
    const CFuint nbNodes = face->nbNodes();
    const vector<Node*>& nodes = *face->getNodes();
    assert(nbNodes <= 4);
    for (CFuint in = 0; in < nbNodes; ++in) {
      mapNodeIDToTrsFaceID.insert(nodes[in]->getGlobalID(), iFace);
    }
    faceBuilder.releaseGE();
  } 
  mapNodeIDToTrsFaceID.sortKeys();
  
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  SafePtr<TopologicalRegionSet> cellTRS = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = states.size();
  assert(nbCells == cellTRS->getLocalNbGeoEnts());
  
  // compute the size of the data needed for describing all partition cells
  // store the corresponding cell local IDs
  vector<CFuint> nbPartitionCellsAndSize(2, 0);
  vector<CFuint> partitionCellsIDs; 
  partitionCellsIDs.reserve(nbCells);
  for (CFuint s = 0; s < nbCells; ++s){
    // only parallel updatable cells can have internal faces that are partition faces for other processors
    // if (states[s]->isParUpdatable()) {
    partitionCellsIDs.push_back(s);
    nbPartitionCellsAndSize[1] += ElementDataArray<0>::getElementSize(cellTRS->getNbNodesInGeo(s), 1);
    // }
  }
  
  nbPartitionCellsAndSize[0] = partitionCellsIDs.size();
  assert(nbPartitionCellsAndSize[0] > 0);
  assert(nbPartitionCellsAndSize[1] > 0);
  assert(nbPartitionCellsAndSize[0] < nbPartitionCellsAndSize[1]);
  
  // compute the maximum sizes for the number of partition cells and their description data
  vector<CFuint> maxNbPartitionCellsAndSize(2, 0);
  
  const std::string nsp = this->getMethodData().getNamespace();
  const CFuint nbPr = PE::GetPE().GetProcessorCount(nsp);
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  const CFuint rank = PE::GetPE().GetRank(nsp);
  MPI_Allreduce(&nbPartitionCellsAndSize[0], &maxNbPartitionCellsAndSize[0], 2, 
                MPIStructDef::getMPIType(&nbPartitionCellsAndSize[0]), MPI_MAX, comm);
  assert(maxNbPartitionCellsAndSize[0] > 0);
  assert(maxNbPartitionCellsAndSize[1] > 0);
  assert(maxNbPartitionCellsAndSize[0] < nbPartitionCellsAndSize[1]);
    
  // store all data needed for the local partition cell description
  ElementDataArray<0> partitionCells;
  partitionCells.reserve(maxNbPartitionCellsAndSize[0], maxNbPartitionCellsAndSize[1]);
  assert(maxNbPartitionCellsAndSize[0] < maxNbPartitionCellsAndSize[1]);
  
  // build only cells belonging to the partition and store all related data
  GeometricEntityPool<CellTrsGeoBuilder> cellBuilder;
  cellBuilder.setup();
  cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellData = cellBuilder.getDataGE();
  cellData.trs = cellTRS;
  
  for (CFuint i = 0; i < partitionCellsIDs.size(); ++i) {
    cellData.idx = partitionCellsIDs[i];
    GeometricEntity *const cell = cellBuilder.buildGE();
    
    // insert partition cell data
    partitionCells.setBeginEptr();
    partitionCells.addElemDataEntry(cell->getState(0)->getGlobalID());
    
    // local ID  (no need to give a consistent value here)
    partitionCells.addElemDataEntry(cell->getState(0)->getLocalID());
    // entity type ID (no need to give a consistent value here)
    partitionCells.addElemDataEntry(0);
    
    const CFuint nbCellNodes = cell->nbNodes();
    partitionCells.addElemDataEntry(nbCellNodes);
    partitionCells.addElemDataEntry(1); // nb states in cell
    for (CFuint n = 0; n < nbCellNodes; ++n) {
      partitionCells.addElemDataEntry(cell->getNode(n)->getGlobalID()); 
    }
    partitionCells.addElemDataEntry(cell->getState(0)->getGlobalID());
    partitionCells.setEndEptr();
    
    cellBuilder.releaseGE();
  }
    
  // array of flags to detect which partition faces have been successfully processed 
  vector<bool> foundID(nbFaces, false);
  
  ElementDataArray<0> commPartitionCells;
  commPartitionCells.resize(maxNbPartitionCellsAndSize[0], maxNbPartitionCellsAndSize[1]);
  
  int ln[3];
  ln[0] = maxNbPartitionCellsAndSize[0];
  ln[1] = maxNbPartitionCellsAndSize[1];
  ln[2] = 1;
  
  for (CFuint root = 0; root < nbPr; ++root) {
    
    CFuint commNbCells = 0;
    // only root copy its data
    if (root == rank) {
      commPartitionCells.copy(partitionCells, false);
      commNbCells = partitionCellsIDs.size();
    }
    
    MPIStruct ms;
    MPIStructDef::buildMPIStruct(commPartitionCells.startPtr(), commPartitionCells.startData(), &commNbCells, ln, ms);
    MPI_Bcast(ms.start, 1, ms.type, root, comm);
    
    // partition face is surely attached to one cell belonging to the partition region and not belonging to the current partition 
    // find such cell by scanning all faces of communicated partition cells
    if (root != rank) {
      ElementDataArray<0>::Itr itr;
      for (CFuint iElem = 0; iElem < commNbCells; ++iElem) {
	itr = commPartitionCells.getElement(iElem); 
	// check if at least one node of the current partition cell belongs to local partition faces 
	bool isFound =  false;
	
	const CFuint nbNodesInCell = itr.get(ElementDataArray<0>::NB_NODES);	
	assert(nbNodesInCell <= 8);	
	for (CFuint jn = 0; jn < nbNodesInCell; ++jn) {
	  pair<MapIt, MapIt> localFaceIDs = mapNodeIDToTrsFaceID.find(itr.getNode(jn), isFound);
	  
	  if (isFound) {
	    // don't exit the loop but process all available faces: the same cell could be connected to multiple partition faces 
	    for (MapIt faceItr = localFaceIDs.first; faceItr != localFaceIDs.second; ++faceItr) {
	      const CFuint faceID = faceItr->second;
	      
	      if (!foundID[faceID]) {
		faceData.idx = faceID;
		GeometricEntity *const face = faceBuilder.buildGE();
		const CFuint nbNodesInFace = face->nbNodes();
		assert(nbNodesInFace < nbNodesInCell);
		
		CFuint counter = 0;
		for (CFuint in = 0; in < nbNodesInFace; ++in) {
		  for (CFuint jn = 0; jn < nbNodesInCell; ++jn) {
		    if (itr.getNode(jn) == face->getNode(in)->getGlobalID()) {
		      counter++;
		      break;
		    }
		  }
		}
		
		// if all global node IDs of the face are contained in the partition cell and the global cell ID is different
		// from the left state of the face, set the global ID in the ghost state of the current partition face
		if (counter == nbNodesInFace && (itr.get(ElementDataArray<0>::GLOBAL_ID) != face->getState(0)->getGlobalID())) {
		  face->getState(1)->setGlobalID(itr.get(ElementDataArray<0>::GLOBAL_ID));
		  foundID[faceID] = true;
		  rankPartitionFaces[faceID] = root;
		}
		
		faceBuilder.releaseGE();
	      } 
	    }
	  }
	}
      }
    }
  }
  
  // sanity check
  for (CFuint i = 0; i < nbFaces; ++i) {
    if (!foundID[i]) {
      cout << "WATCH OUT: global ghost state ID not found for Face [" << i << "] in P-" << rank  << endl;
    }
  }
  
#endif
}

//////////////////////////////////////////////////////////////////////////////

void StencilCUDASetup::computeStencil()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  // table that stores the connectivity state-neighbor states
  vector< vector<State*> > neighborStates(nbStates);
  
  SelfRegistPtr<ComputeStencil> computeStencil(Environment::Factory<ComputeStencil>::getInstance().
					       getProvider(_stencilType)->create(_stencilType));
  configureNested ( computeStencil.getPtr(), m_stored_args );
  
  computeStencil->setDataSocketSinks(socket_states, socket_nodes, socket_stencil, socket_gstates);
  (*computeStencil)();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
