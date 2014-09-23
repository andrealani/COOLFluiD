#include <memory>

#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"
#include "Framework/RadiationLibrary.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"

#include "RadiativeTransferSanna/RadiativeTransferSanna.hh"
#include "RadiativeTransferSanna/Slab1DFVMCC.hh"
#include "FiniteVolume/CellCenterFVM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Slab1DFVMCC, DataProcessingData, RadiativeTransfer> 
slab1DFVMCCProvider("Slab1DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >("RadiationLibrary","Name of the radiation library."); 
  options.addConfigOption< vector<CFreal> >("StagnationPoint","Coordinates (x,y,z) of stagnation point.");
}
      
//////////////////////////////////////////////////////////////////////////////

Slab1DFVMCC::Slab1DFVMCC(const std::string& name) :
  DataProcessingCom(name),
  m_radLibrary(),
  socket_qrad("qrad"),
  socket_qradFluxWall("qradFluxWall"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_gstates("gstates"),
  socket_normals("normals"),
  m_stagPoint(),
  m_nbWallFaces(0),
  m_nodeIdToStateId(),
  m_orderedFaceIDs(),
  m_stagnationLineCells(CFNULL),
  m_fvmccData(CFNULL),
  m_meshByLines(),
  m_deltal()
{
  addConfigOptionsTo(this);
 
  m_radLibraryName = "Parade";
  setParameter("RadiationLibrary",&m_radLibraryName);
  
  m_stagnationPointXYZ = vector<CFreal>();
  setParameter("StagnationPoint",&m_stagnationPointXYZ);
}

//////////////////////////////////////////////////////////////////////////////

Slab1DFVMCC::~Slab1DFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCC::unsetup()
{
  for (CFuint i = 0; i < m_meshByLines.size(); ++i) {
    deletePtr(m_meshByLines[i]);
  }
  
  for (CFuint i = 0; i < m_deltal.size(); ++i) {
    deletePtr(m_deltal[i]);
  }
}
 
//////////////////////////////////////////////////////////////////////////////
     
std::vector<Common::SafePtr<BaseDataSocketSink> > Slab1DFVMCC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);
  
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > Slab1DFVMCC::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_qrad);
  result.push_back(&socket_qradFluxWall);
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////
      
void Slab1DFVMCC::setup()
{
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  m_nodeIdToStateId.resize(states.size());
  for (CFuint i = 0; i < m_nodeIdToStateId.size(); ++i) {
    m_nodeIdToStateId[i] = i;
  }
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // default stagnation point is NOT set
  if (m_stagnationPointXYZ.size() > 0) {
    cf_always_assert(m_stagnationPointXYZ.size() == dim);
   
    m_stagPoint.resize(dim);
    for (CFuint i = 0; i < dim; ++i) {
      m_stagPoint[i] = m_stagnationPointXYZ[i];
    }
  }
  else {
    m_stagPoint.resize(dim);
    m_stagPoint = 0.0;
  }
  
  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());
  
  m_fvmccData = fvmcc->getData();
  
  //  m_nsVarSet = m_fvmccData->getDiffusiveVar().d_castTo<NavierStokes2DVarSet>();
  
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts(); //
  
  DataHandle<CFreal> qrad = socket_qrad.getDataHandle(); //getGlobalNbGeoEnts()
  qrad.resize(nbCells);
  
  // build lines of cells starting from the wall and marching towards the outer boundaries
  
  // face builder 
  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > 
    faceBuilder = m_fvmccData->getFaceTrsGeoBuilder();
  
  SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
  faceBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  FaceTrsGeoBuilder::GeoData& faceData = faceBuilder->getDataGE();
  faceData.isBFace = true;
  
  // cell builder
  SafePtr<GeometricEntityPool<CellTrsGeoBuilder> > cellBuilder = m_fvmccData->getCellTrsGeoBuilder();
  
  SafePtr<CellTrsGeoBuilder> cellBuilderPtr = cellBuilder->getGeoBuilder();
  cellBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  CellTrsGeoBuilder::GeoData& cellData = cellBuilder->getDataGE();
  cellData.trs = cells;
    
  m_nbWallFaces = 0;
  
  vector< Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
  
  if (trsList.size() != 1) {cout << "Only one Wall boundary TRS is supported for now"<< endl; abort();}
  
  for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    m_nbWallFaces += trsList[iTRS]->getLocalNbGeoEnts();
  } 
  
  // allocate the memory for the radiative heat flux at the wall 
  socket_qradFluxWall.getDataHandle().resize(m_nbWallFaces);
  
  // allocate the full memory for the mapping faceIDs to lines of cells
  m_meshByLines.reserve(m_nbWallFaces);
  
  // dx in line direction
  m_deltal.reserve(m_nbWallFaces);
  
  m_orderedFaceIDs.reserve(m_nbWallFaces);
  
  CFMap<CFreal, CFuint> mapXToFaceID;
  mapXToFaceID.reserve(m_nbWallFaces);
  
  // start from the inner boundary and march towards the outer boundary
  // outer boundary is detected when correspodning face has one ghost state
  CFuint maxNbCellsInLine = 0; 
  CFreal minDistanceFromStagPoint = 10000000000.; // starting huge value
  
  RealVector faceCenter0(dim);
  RealVector faceCenter1(dim);
  
  for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trsList[iTRS];
    faceData.trs = currTrs;
    
    const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
      // build the GeometricEntity
      faceData.idx = iFace;
      GeometricEntity* const currFace = faceBuilder->buildGE();
      const CFuint startFaceID = currFace->getID();
      
      CFuint faceID = startFaceID;
      CFuint cellID = currFace->getState(LEFT)->getLocalID();
      bool outerBoundaryReached = false;
      
      vector<CFuint>* cellLines = new vector<CFuint>();
      vector<CFreal>* deltal = new vector<CFreal>();
      // preallocate the memory for each line if you already computed for previous faces
      if (maxNbCellsInLine > 0) {
	cellLines->reserve(maxNbCellsInLine);
	deltal->reserve(maxNbCellsInLine);
      }
      
      // assume 2D 
      faceCenter0 = 0.5*(*currFace->getNode(0) + *currFace->getNode(1));
      
      while (!outerBoundaryReached) {
	// build cell corresponding to the innerState
	cellData.idx = cellID;
	
	// store the current cell ID
	cellLines->push_back(cellID);
	
	GeometricEntity* const cell = cellBuilderPtr->buildGE();
	
	// keep track of the minimum distance to stagnation point and corresponding faceID
	if (computeOnStagnationLine()) {
	  const CFreal distance = MathFunctions::getDistance(cell->getState(0)->getCoordinates(), m_stagPoint);
	  if (distance < minDistanceFromStagPoint) {
	    minDistanceFromStagPoint = distance;
	    // set the pointer to the stagnation line cells
	    m_stagnationLineCells = cellLines;
	  }
	}
	
	const vector<GeometricEntity*>& cellFaces = *cell->getNeighborGeos(); 
	const CFuint nbCellFaces = cellFaces.size(); 
	
	// determine the local idx (0-nbCellFaces) of the face in the current cell
	CFuint cellFaceIdx = 0;
	bool cellFaceIdxFound = false;
	for (CFuint f = 0; f < nbCellFaces; ++f) {
	  if (cellFaces[f]->getID() == faceID) {
	    cellFaceIdx = f;
	    cellFaceIdxFound = true;
	    break;
	  }
	}
	cf_assert(cellFaceIdxFound);
	
	// find the next cell ID as the cell having ID != current cell ID, being on one side of the opposite face ID
	const CFuint oppositeCellFaceIdx = m_fvmccData->getOppositeIFace(cellFaceIdx, dim, cell->nbNodes());
	GeometricEntity *const oppositeFace = cellFaces[oppositeCellFaceIdx];
	faceCenter1 = 0.5*(*oppositeFace->getNode(0) + *oppositeFace->getNode(1));
	
	const CFreal distanceF0F1 = MathFunctions::getDistance(faceCenter0, faceCenter1);
	cf_assert(distanceF0F1 > 0.);
	deltal->push_back(distanceF0F1);
	
	// face center 1 becomes the new face center 0 
	faceCenter0 = faceCenter1;
	
	  // terminate the algorithm if you reach another boundary (outer) face 
	if (oppositeFace->getState(1)->isGhost()) {
	  maxNbCellsInLine = std::max(maxNbCellsInLine, static_cast<CFuint>(cellLines->size()));
	  
	  // insert the full storage of cell lines starting from the given face
	  m_meshByLines.push_back(cellLines);
	  m_deltal.push_back(deltal);
	  cf_assert(cellLines->size() == deltal->size());
	  
	  // compute mid point of face
	  const CFuint nbFaceNodes = oppositeFace->nbNodes();
	  CFreal xcoord = 0.0;
	  for (CFuint n = 0; n < nbFaceNodes; ++n) {
	    xcoord += (*oppositeFace->getNode(n))[XX];
	  }
	  xcoord /= nbFaceNodes;
	  
	  // insert the pair <x coordinate, faceID>
	  mapXToFaceID.insert(xcoord, m_meshByLines.size() -1);
	  outerBoundaryReached = true;
	}
	else {
	  const CFuint nextCellID = (oppositeFace->getState(0)->getLocalID() == cellID) ?  
	    oppositeFace->getState(1)->getLocalID() : oppositeFace->getState(0)->getLocalID();
	  cf_assert(oppositeFace->getState(0)->getLocalID() == cellID || 
		    oppositeFace->getState(1)->getLocalID() == cellID);
	  
	  //assign the nextCellID to the current cellID
	  cellID = nextCellID;
	  // assign the next face ID to the current face ID and restart the whole process
	  faceID = cellFaces[oppositeCellFaceIdx]->getID();
	}
	
	cellBuilderPtr->releaseGE();
      }
      
      // release the geometric entity
      faceBuilder->releaseGE();
    }  
  }
  
  // sort the entries in the map
  mapXToFaceID.sortKeys();
  
  cf_always_assert(mapXToFaceID.size() == m_nbWallFaces);
  
  // we keep track of the ordering of faces [0, m_nbWallFaces - 1] based on (growing) X coordinate
  for (CFuint i = 0; i < m_nbWallFaces; ++i) {
    m_orderedFaceIDs[i] = mapXToFaceID[i];
  }
  
  // reorder the mesh by lines following the X direction
  vector<vector<CFuint>* > bkp(m_nbWallFaces);
  bkp = m_meshByLines;
  for (CFuint i = 0; i < m_nbWallFaces; ++i) {
    m_meshByLines[i] = bkp[m_orderedFaceIDs[i]];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);  
  
  m_radLibrary = Environment::Factory<RadiationLibrary>::getInstance().
    getProvider(m_radLibraryName)->create(m_radLibraryName);
  
  cf_assert(m_radLibrary.isNotNull());
  
  configureNested ( m_radLibrary.getPtr(), args );
}
      
//////////////////////////////////////////////////////////////////////////////

void Slab1DFVMCC::executeOnTrs()
{
  // careful here !!! one TRS at a time !!!!
  cf_assert(getTrsList().size() == 1);
  
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  auto_ptr<ProxyDofIterator<CFreal> > nstatesProxy
    (new DofDataHandleIterator<CFreal, State, GLOBAL>(states,&m_nodeIdToStateId));
  
  if (computeOnStagnationLine()) {
    m_radLibrary->runOnStagnationLine(m_stagnationLineCells,
				      nstatesProxy.get(),
				      &socket_qrad.getDataHandle()[0]); 
  }
  else {
    m_radLibrary->runOnStructuredMesh(m_meshByLines,
				      nstatesProxy.get(),
				      &socket_qrad.getDataHandle()[0]); 
    
    computeWallRadHeatFlux();
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void Slab1DFVMCC::computeWallRadHeatFlux()
{
  DataHandle<CFreal> qrad = socket_qrad.getDataHandle();
  DataHandle<CFreal> qradFluxWall = socket_qradFluxWall.getDataHandle();
  
  for (CFuint iFace = 0; iFace < m_meshByLines.size(); ++iFace) {
    const CFint end = m_meshByLines[iFace]->size() - 1;
    for (CFint il = end; il >= 0; il--) {
      const CFuint trsFaceID = m_orderedFaceIDs[iFace];
      const CFuint cellID    = (*m_meshByLines[iFace])[il];
      qradFluxWall[trsFaceID] += qrad[cellID]*(*m_deltal[trsFaceID])[il];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
