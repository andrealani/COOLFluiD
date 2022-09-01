#include "Common/PE.hh"
#include "Common/EventHandler.hh"
#include "MathTools/LeastSquaresSolver.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

#include "Framework/CFSide.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/TRSDistributeData.hh"
#include "Framework/State.hh"

#include "Framework/LSSData.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <iostream>
#include <limits>
#include <vector>
#include <typeinfo>

#include "FluxReconstructionMethod/FluxReconstructionSolver.hh" 

#include "FluxReconstructionMethod/FluxReconstruction.hh"

#include "FluxReconstructionMethod/MeshFittingAlgorithmFRQ2Quads.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"


#include "FiniteVolume/CellData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"
#include "MeshTools/ComputeWallDistanceVectorFRMPI.hh"
#include "MeshTools/MeshToolsFR.hh"



#include "Common/CFMultiMap.hh"
#include <math.h>
#include <cmath>
#include <fstream>

#include "Common/PE.hh"

#ifdef CF_HAVE_MPI
#include "Common/MPI/MPIStructDef.hh"
#endif


///////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

///////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
        
///////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshFittingAlgorithmFRQ2Quads,
		      DataProcessingData,
		      FluxReconstructionModule>
MeshFittingAlgorithmFRQ2QuadsFluxReconstructionProvider("MeshFittingAlgorithmFRQ2Quads");

///////////////////////////////////////////////////////////////////


void MeshFittingAlgorithmFRQ2Quads::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("minPercentile","Percentile for minimum spring value");
  options.addConfigOption< CFreal >("maxPercentile","Percentile for maximum spring value");
  options.addConfigOption< CFreal >("meshAcceleration","How fast the mesh moves in mesh steps");
  options.addConfigOption< CFuint >("monitorVarID","Monitor variable ID (from State) for mesh adaptation");
  options.addConfigOption< CFuint >("monitorPhysVarID","Monitor physical variable ID (from physical data) for mesh adaptation");
  options.addConfigOption< CFreal >("equilibriumSpringLength","Length of spring for equilibrium");
  options.addConfigOption< CFreal >("ratioBoundaryToInnerEquilibriumSpringLength","ratio between the equilibrium length of a Boundary spring to an Inner spring");
  options.addConfigOption< std::vector<std::string> >("unlockedBoundaryTRSs","TRS's to be unlocked");
  options.addConfigOption< CFreal >("AcceptableDistance","Distance from user-defined boundary");
  options.addConfigOption< bool   >("ThetaMid","Semi torsional Sping analogy for 2D quadrilateral mesh based on the middle egde-facing  angle (true) or the 3 edge-facing angles (false) ");
  options.addConfigOption< bool   >("InterpolateState"," State Interplation to dissociate the nodal movement and the solution in each CC. ");


  options.addConfigOption< bool >("smoothSpringNetwork","smooth the spring network");
  options.addConfigOption< bool >("smoothNodalDisp","smooth the nodal displacememt");
}

///////////////////////////////////////////////////////////////////

MeshFittingAlgorithmFRQ2Quads::MeshFittingAlgorithmFRQ2Quads(const std::string& name) :

  Framework::DataProcessingCom(name),  
  socket_stiffness("stiffness"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_normals("normals"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  socket_rhs("rhs"), 
  socket_wallDistance("wallDistance",false),
  socket_nodeisAD("nodeisAD", false),
  socket_nodeDistance("nodeDistance", false),
  //socket_stencil("stencil"),
  m_lss(CFNULL),
  m_geoBuilder(),
  m_cellBuilder(CFNULL),
  m_faceBuilder(),
  m_order(),
  m_frData(CFNULL),
  m_solPolyValsAtNodes(CFNULL),
  m_nbrNodesElem(),
  m_nbrSolPnts(),
  m_nbrFaceFlxPnts(),
  m_face(),
  m_cells(),
  m_states(),
  m_faceFlxPntConnPerOrient(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_faceConnPerOrient(CFNULL),
  m_orient(),
  m_cellStatesFlxPnt(),
  m_nbrSolDep(),
  m_flxSolDep(),
  m_solPolyValsAtFlxPnts(),
  m_flxLocalCoords(CFNULL),
  m_flxPntCoords(),
  m_nbrEqs(),
  m_dim(),
  m_flxPntsLocalCoords(),
  m_unitNormalFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_allCellFlxPnts(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_faceJacobVecs(),
  m_normalsAMR(),
  m_orientsAMR(),
  m_intCell(),
  m_cellStates(),
  m_cellStatesFlxPntBnd(),
  #ifdef CF_HAVE_MPI
  m_comm(MPI_COMM_WORLD),
  #endif
  m_myRank(0),
  m_nbProc(1),
  m_iElemType(),
  m_cell(),
  oldStates(),
  oldCoordinates(),
  m_firstLeftState(),
  m_secondLeftState(),
  m_firstRightState(),
  m_secondRightState()
{
  this->addConfigOptionsTo(this);
  
  m_minPercentile = 0.10;
  this->setParameter("minPercentile", &m_minPercentile);

  m_maxPercentile = 0.90;
  this->setParameter("maxPercentile", &m_maxPercentile);

  m_meshAcceleration = 0.1;  
  this->setParameter("meshAcceleration", &m_meshAcceleration);

  m_monitorVarID = std::numeric_limits<CFuint>::max();
  this->setParameter("monitorVarID", &m_monitorVarID);
  
  m_monitorPhysVarID = std::numeric_limits<CFuint>::max();
  this->setParameter("monitorPhysVarID", &m_monitorPhysVarID);
  
  m_equilibriumSpringLength = 0.;
  this->setParameter("equilibriumSpringLength", &m_equilibriumSpringLength);

  m_ratioBoundaryToInnerEquilibriumSpringLength = 1.;
  this->setParameter("ratioBoundaryToInnerEquilibriumSpringLength", &m_ratioBoundaryToInnerEquilibriumSpringLength);
  
  this->setParameter("unlockedBoundaryTRSs", &m_unlockedBoundaryTRSs);

  m_acceptableDistance = 0.;
  this->setParameter("AcceptableDistance",&m_acceptableDistance);

  m_thetaMid = true;
  this->setParameter("ThetaMid",&m_thetaMid);

  m_interpolateState = false;
  this->setParameter("InterpolateState",&m_interpolateState);

  m_smoothSpringNetwork = false;
  this->setParameter("smoothSpringNetwork",&m_smoothSpringNetwork);

  m_smoothNodalDisp = false;
  this->setParameter("smoothNodalDisp",&m_smoothNodalDisp);

  
  
}

//////////////////////////////////////////////////////////////////////////////

MeshFittingAlgorithmFRQ2Quads::~MeshFittingAlgorithmFRQ2Quads()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

 Framework::DataProcessingCom::configure(args);
  
  cf_assert(m_minPercentile >= 0.);
  cf_assert(m_minPercentile < m_maxPercentile);
  cf_assert(m_maxPercentile <= 1.);
  cf_assert(m_meshAcceleration > 0. && m_meshAcceleration < 1.);
  cf_assert(m_acceptableDistance >= 0);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
MeshFittingAlgorithmFRQ2Quads::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_normals);
  result.push_back(&socket_rhs);
  result.push_back(&socket_wallDistance);
  result.push_back(&socket_nodeisAD);  
  result.push_back(&socket_nodeDistance);
  //result.push_back(&socket_stencil);
  return result;
}
//////////////////////////////////////////////////////////////////////////////
std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
MeshFittingAlgorithmFRQ2Quads::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_stiffness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::setup()
{
  CFAUTOTRACE;
  
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  DataProcessingCom::setup();

   // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  //resize and initialize the storage of the nodal stiffness
  DataHandle< CFreal > stiffness = socket_stiffness.getDataHandle();
  stiffness.resize(socket_nodes.getDataHandle().size());
  stiffness = 0.;

    // AL: this might be useless ... (@see MeshRigidMove/StdSetup.cxx)
  SubSystemStatusStack::getActive()->setMovingMesh(true);
   
  //m_lss = getMethodData().getLinearSystemSolver()[0];
  const std::string name = getMethodData().getNamespace();
  m_lss = getMethodData().getCollaborator<LinearSystemSolver>(name);  

  //CFLog(VERBOSE, "MeshFittingAlgorithmFRQ2Quads::setup() -----=> LSS is " << m_lss->getName() << "\n");
  
  Common::SafePtr<Framework::SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();

  Common::SafePtr<FluxReconstructionSolver> frsolver = spaceMethod.d_castTo<FluxReconstructionSolver>(); //# Change here
  cf_assert(frsolver.isNotNull());
  m_frData = frsolver->getData();
  

  vector< FluxReconstructionElementData* >& frLocalData = m_frData->getFRLocalData();

 
  getMethodData().getUpdateVarSet()->setup();


  if (m_monitorVarID == std::numeric_limits<CFuint>::max() ||
    m_monitorVarID > PhysicalModelStack::getActive()->getNbEq()) {
    CFLog(WARN, "MeshFittingAlgorithmFRQ2Quads::setup() => monitorVarID not specified or invalid: will be set to 0\n");
    m_monitorVarID = 0;
  }
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);

  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();

  m_order = static_cast<CFuint>(order);


  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtNodes = frLocalData[0]->getCoefSolPolyInNodes();

  m_nbrNodesElem = m_solPolyValsAtNodes->size();
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  m_cellBuilder = m_frData->getCellBuilder();
  m_geoBuilder =  m_frData->getStdTrsGeoBuilder();
  m_faceBuilder = m_frData->getFaceBuilder();

  m_cells.resize(2);
  m_states.resize(2);
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = flxLocalCoords->size();
  
  // number of sol points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  
  cf_assert(m_nbrSolPnts == (frLocalData[0]->getSolPntsLocalCoords())->size());

  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

   // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();


    // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  m_flxPntCoords.resize(m_nbrFaceFlxPnts);

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
  }
  
  // create internal and ghost states
  m_cellStatesFlxPnt.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt[LEFT].push_back(new State());
    m_cellStatesFlxPnt[RIGHT].push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt[LEFT][iFlx]->setLocalID(iFlx);
    m_cellStatesFlxPnt[RIGHT][iFlx]->setLocalID(iFlx);
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPntBnd.push_back(new State());
    m_cellStatesFlxPntBnd.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPntBnd[iFlx]->setLocalID(iFlx);
    m_cellStatesFlxPntBnd[iFlx]->setLocalID(iFlx);
  }

  

  m_flxSolDep = frLocalData[0]->getFlxPntSolDependency();

  m_nbrSolDep = ((*m_flxSolDep)[0]).size();

  m_flxPntsLocalCoords.resize(m_nbrFaceFlxPnts);
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecs.resize(m_nbrFaceFlxPnts);

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_faceJacobVecs[iFlx].resize(m_dim);
  }
  // get all flux points of a cell
  m_allCellFlxPnts = frLocalData[0]->getFlxPntsLocalCoords();

  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDir();

  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
  
  createGeneralConnectivityFR(); 

  createNodalConnectivity();

  CFuint nbTotalFaces = 0;

  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");
  nbTotalFaces = faces->getLocalNbGeoEnts();

   #ifdef CF_HAVE_MPI
  const std::string nsp = this->getMethodData().getNamespace();
  m_comm   = PE::GetPE().GetCommunicator(nsp);
  m_myRank = PE::GetPE().GetRank(nsp);
  m_nbProc = PE::GetPE().GetProcessorCount(nsp);
  #endif

  // loop over all processors
  for (CFuint root = 0; root < m_nbProc; ++root) {
	for(CFuint iTrs=0; iTrs<m_unlockedBoundaryTRSs.size();++iTrs){
    		FaceToCellGEBuilder::GeoData& facesData = m_faceBuilder->getDataGE();
    		m_faceBuilder->getDataGE().isBoundary = true;
    		Common::SafePtr<Framework::TopologicalRegionSet> wallFaces =
   	         Framework::MeshDataStack::getActive()->getTrs( m_unlockedBoundaryTRSs[iTrs] );
    		facesData.facesTRS = wallFaces;
    		nbTotalFaces += wallFaces->getLocalNbGeoEnts();
    		map< std::string , vector< vector< CFuint > > >&
    		bndFacesStartIdxsPerTRS = m_frData->getBndFacesStartIdxs();
    		vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[wallFaces->getName()];
     		cf_assert(bndFacesStartIdxs.size() != 0);
		// number of face orientations (should be the same for all TRs)
     		const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

     		// number of TRs
     		const CFuint nbTRs = wallFaces->getNbTRs();
     		cf_assert(bndFacesStartIdxs.size() == nbTRs);
  	}
  }

  m_normalsAMR.resize(nbTotalFaces*m_dim);
  m_orientsAMR.resize(nbTotalFaces*m_dim);
  
  
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();

  findBoundaryNodes();

  for (CFuint root = 0; root < m_nbProc; ++root) {
  	computeMovingInBoundaryNodeNormals();  
  }

  DataHandle<bool> nodeisAD = socket_nodeisAD.getDataHandle();
  nodeisAD.resize(socket_nodes.getDataHandle().size());
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
 	nodeisAD = false;
  }

  oldStates.resize(m_nbrEqs*socket_states.getDataHandle().size());
  oldStates = 0.;

  oldCoordinates.resize(m_dim*socket_states.getDataHandle().size());
  oldCoordinates = 0.;

  m_firstLeftState.resize(m_nbrEqs);
  m_secondLeftState.resize(m_nbrEqs);
  m_firstRightState.resize(m_nbrEqs);
  m_secondRightState.resize(m_nbrEqs);

  if(m_interpolateState){
    CFLog(INFO, "Interpolation Activated ---- Save Old Mesh properties\n");
    saveOldMeshProperties();
  }
}

//////////////////////////////////////////////////////////////////////////////


void MeshFittingAlgorithmFRQ2Quads::unsetup()
{
  CFAUTOTRACE;
  Framework::DataProcessingCom::unsetup();
}
  

//////////////////////////////////////////////////////////////////////////////


 void MeshFittingAlgorithmFRQ2Quads::createNodalConnectivity()
{ 
  /////////////////////////
  //2D quadrilateral mesh//
  /////////////////////////

   // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
  Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  
  if( nbElemTypes==4 && nbDims==2){

	  // Connectivity information 2D quadrilateral
	  CFuint nbPairsNodeNode = 200000; 
	  
	  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
	  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeNode(nbPairsNodeNode);
	  std::multimap<CFuint, CFuint>  mapNodeNode;

	  typedef CFMultiMap<CFuint, CFuint > MapFaceNode;
	  Common::CFMultiMap<CFuint,CFuint>  m_mapFaceNode(200000);
	  typedef MapFaceNode::MapIterator mapFaceIt;

	  typedef CFMultiMap<CFuint, CFuint > MapNodeFace;
	  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeFace(200000);
	  typedef MapNodeFace::MapIterator mapNodeIt;

	  typedef CFMultiMap<CFuint, CFuint > MapNodeTRS;
	  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeTRS(200000);
	  typedef MapNodeTRS::MapIterator mapNodeTRSIt;

	  typedef CFMultiMap<CFuint, CFuint > MapNodeCell;
	  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeCell(200000);
	  typedef MapNodeCell::MapIterator mapNodeCellIt;
	  
	  const CFuint nbCells = cells->getLocalNbGeoEnts();
	  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
	  geoData.trs = cells;
	  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");
	  
	  // get the face start indexes
	  vector< CFuint >& innerFacesStartIdxs = m_frData->getInnerFacesStartIdxs();
	  // get number of face orientations
	  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

	  // get the geodata of the face builder and set the TRSs
	  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
	  geoDataFace.cellsTRS = cells;
	  geoDataFace.facesTRS = faces;
	  geoDataFace.isBoundary = false;

	  
	  // loop over different orientations
	  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
	  {
	    // start and stop index of the faces with this orientation
	    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
	    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

	    // loop over faces with this orientation
	    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
	    {
	      // build the face GeometricEntity
	      geoDataFace.idx = faceID;
	      m_face = m_faceBuilder->buildGE();
	      std::vector<Framework::Node*>& faceNodes = *m_face->getNodes();
	      for (CFuint iNodeC = 0; iNodeC < faceNodes.size(); ++iNodeC){
			m_mapFaceNode.insert(faceID,faceNodes[iNodeC]->getLocalID());
			m_mapNodeFace.insert(faceNodes[iNodeC]->getLocalID(), faceID);
	      }
	      m_mapNodeNode.insert(faceNodes[0]->getLocalID(),faceNodes[1]->getLocalID());
	      m_mapNodeNode.insert(faceNodes[1]->getLocalID(),faceNodes[0]->getLocalID());
	      m_mapFaceNode.sortKeys();
	      m_mapNodeFace.sortKeys();
	      m_mapNodeNode.sortKeys();
	      m_faceBuilder->releaseGE();
	    }
	   }

	  const CFuint nbrElemTypes = elemType->size();
	  cf_assert(nbrElemTypes == 1);
	  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
	  {
	    // get start and end indexes for this type of element
	    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
	    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

	    // loop over cells
	    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
	    {
		geoData.idx = elemIdx;
		m_cell = m_cellBuilder->buildGE();
		std::vector<Framework::Node*>& cellNodes = *m_cell->getNodes();
		for (CFuint iNodeC = 0; iNodeC < cellNodes.size(); ++iNodeC) // loop over nodes
		{
			m_mapNodeCell.insert(cellNodes[iNodeC]->getLocalID(),elemIdx);
			
		} 
	   	m_mapNodeCell.sortKeys();
		m_cellBuilder->releaseGE();
	    }
	  }
	  
	  m_mapNodeCell1 = m_mapNodeCell;
	  
	   //loop over unlocked boundaries
	  for(CFuint iTrs=0; iTrs<m_unlockedBoundaryTRSs.size();++iTrs){
	    FaceToCellGEBuilder::GeoData& facesData = m_faceBuilder->getDataGE();
	    m_faceBuilder->getDataGE().isBoundary = true;
	    Common::SafePtr<Framework::TopologicalRegionSet> wallFaces =
	    Framework::MeshDataStack::getActive()->getTrs( m_unlockedBoundaryTRSs[iTrs] );
	    facesData.facesTRS = wallFaces;
	    map< std::string , vector< vector< CFuint > > >&
	    bndFacesStartIdxsPerTRS = m_frData->getBndFacesStartIdxs();
	    vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[wallFaces->getName()];
	    // number of face orientations (should be the same for all TRs)
	    cf_assert(bndFacesStartIdxs.size() != 0);
	    const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;
	    // number of TRs
	    const CFuint nbTRs = wallFaces->getNbTRs();
	    cf_assert(bndFacesStartIdxs.size() == nbTRs);
	    FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
	    geoData.facesTRS = wallFaces;
	    geoData.isBoundary = true;
	    for (CFuint iTR = 0; iTR < nbTRs; ++iTR){
		  for (m_orient = 0; m_orient < nbOrients; ++m_orient)
		   {
		      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
		      // select the correct flx pnts on the face out of all cell flx pnts for the current orient
		      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
		      {
		        m_flxPntsLocalCoords[iFlx] = (*m_allCellFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
		      }
		      // start and stop index of the faces with this orientation
		      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
		      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

		      // loop over faces with this orientation
		      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
		      {
		        // build the face GeometricEntity
		        geoData.idx = faceID;
		        m_face = m_faceBuilder->buildGE();
		        std::vector<Framework::Node*>& faceNodes = *m_face->getNodes();
		        for (CFuint iNodeC = 0; iNodeC < faceNodes.size(); ++iNodeC){
				m_mapFaceNode.insert(faceID,faceNodes[iNodeC]->getLocalID());
				m_mapNodeFace.insert(faceNodes[iNodeC]->getLocalID(), faceID);
		        }
		        m_mapNodeNode.insert(faceNodes[0]->getLocalID(),faceNodes[1]->getLocalID());
		        m_mapNodeNode.insert(faceNodes[1]->getLocalID(),faceNodes[0]->getLocalID());
		        m_mapNodeTRS.insert(faceNodes[0]->getLocalID(),iTrs);
		        m_mapNodeTRS.insert(faceNodes[1]->getLocalID(),iTrs);
		        m_mapFaceNode.sortKeys();
		        m_mapNodeFace.sortKeys();
		        m_mapNodeNode.sortKeys();
		        m_mapNodeTRS.sortKeys();
		        m_faceBuilder->releaseGE();
		      }
		  }
	     }
	  }
 	m_mapNodeFace1 = m_mapNodeFace;
 	m_mapFaceNode1 = m_mapFaceNode;
 	m_mapNodeNode1 = m_mapNodeNode;
 	m_mapNodeTRS1 = m_mapNodeTRS;
  }
  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR();
  
}
//////////////////////////////////////////////////////////////////////////////*/




//////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::execute()
{
  CFAUTOTRACE;
  CFLog(INFO, "MeshFittingAlgorithmFRQ2Quads::execute() => start \n");
  
  determineIsNodeAD(); // fills the socket nodeisAD depending on the acceptable distance

  resizeSystemSolverToNodalData();
  
  computeNodeStates(); // extrapolates the polynomials at the solution points

  computeSpringTruncationData();

  solveLinearSystem();
  
  updateNodePositions(); // mesh node repositionning
  
  resizeSystemSolverToStateData();
  
  triggerRecomputeMeshData();
  
  CFLog(INFO, "MeshFittingAlgorithm::execute() => end \n");
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::determineIsNodeAD(){
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<bool> nodeisAD = socket_nodeisAD.getDataHandle();
  
  
   for (CFuint iNode = 0; iNode < nodes.size(); ++iNode){
	nodeisAD[nodes[iNode]->getLocalID()] = false;
	if (nodeDistance[nodes[iNode]->getLocalID()] < m_acceptableDistance){
		nodeisAD[nodes[iNode]->getLocalID()] = true;
   	}
  
   }
}

///////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::computeSpringTruncationData()
{
  //This method computes the average, min and max spring constants through a user defined quantile
  //The P^2 algorithm was used to compute the quantiles using the BOOST statistical accumulators
  //This method has very low memory requirements compared to a naive implementation but lower accuracy. 
  //The sum of accumulators are not defined so a parallel reduce can't be used
  //Instead, all the spring constants need to be gathered on process 0 to be accumulated there
  //In particular, the P^2 algorithm will be dependent on the ordering of the springs and thus,
  // on the number of processors and size of the send buffers

  CFAUTOTRACE;
  CFLogDebugMin( "MeshFittingAlgorithm::computeSprings()" << "\n");
  
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle(); 

  using namespace boost::accumulators;
  typedef accumulator_set<CFreal, stats<tag::p_square_quantile> > accumulator_t;
  accumulator_t minQuantileAcc(quantile_probability = m_minPercentile);
  accumulator_t maxQuantileAcc(quantile_probability = m_maxPercentile);
  accumulator_set<double, stats<tag::mean> > meanAcc;

  const std::string nsp = this->getMethodData().getNamespace();
  const int nbProcesses = Common::PE::GetPE().GetProcessorCount(nsp);
  const int processRank = Common::PE::GetPE().GetRank(nsp);
  MPI_Comm communicator = Common::PE::GetPE().GetCommunicator(nsp); 
  
  //Limit the size of the receiving buffer on processor 0
  const CFuint sizeSendBuffer = std::ceil(1000./static_cast<CFreal>(nbProcesses)) + 1;
  const CFuint sizeRecvBuffer = sizeSendBuffer*nbProcesses;

  std::vector<CFreal> sendBuffer, recvBuffer; 
  sendBuffer.reserve(sizeSendBuffer);
  if(processRank == 0) recvBuffer.reserve(sizeRecvBuffer);
  
  SimpleEdgeGraphFR::iterator it = m_edgeGraph.begin();
  bool allSpringsCalculated = false;
  do{
    //compute and store the spring constants in a temporary buffer
    for(CFuint iBuffer=0; iBuffer<sizeSendBuffer && it != m_edgeGraph.end(); ++it){
      Framework::Node* node1 = (*it).firstNode;
     Framework::Node* node2 = (*it).secondNode;
     if(!isNodeLocked(node1) && node1->isParUpdatable()){
     	typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
     	typedef MapNodeNode::MapIterator mapIt;
     	bool found = false;
     	std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node1->getLocalID(), found);
     	cf_assert(found);    
     	for (mapIt it = ite.first; it != ite.second; ++it) {
            	if (node2->getLocalID() == it->second){
        			const CFreal springConstant = computeSpringConstant(node1, node2);
        			if(springConstant > 1.e-4){
          			sendBuffer.push_back(springConstant);
          			++iBuffer;

			}
        		}
      	 }
      }
    }
    //get global number of springs constants in all the buffers
    int sendBufferSize = sendBuffer.size();
    std::vector<int> recvBufferCounts(nbProcesses,0); 
    std::vector<int> recvBufferDisp(nbProcesses,0);
    MPI_Allgather(&sendBufferSize, 1, MPI_INT, &recvBufferCounts[0], 1, MPI_INT, communicator);
    //compute the buffer counts and displacements 
    int sumRecvBufferSize = 0;
    for(CFuint i=0; i< nbProcesses; ++i) {
      recvBufferDisp[i] = sumRecvBufferSize;
      sumRecvBufferSize += recvBufferCounts[i];
    }
    if(processRank == 0) recvBuffer.resize(sumRecvBufferSize); 
    //gather the temporary buffers in processor 0
    //to perform the statistical computations
    MPI_Gatherv(&sendBuffer[0], sendBufferSize, MPI_DOUBLE, &recvBuffer[0], 
        &recvBufferCounts[0], &recvBufferDisp[0], MPI_DOUBLE, 0, communicator); 
    sendBuffer.clear();  
    if(processRank == 0) {
	for(CFuint i=0; i<recvBuffer.size(); ++i){
	const CFreal ke = recvBuffer[i];
         minQuantileAcc(ke);
         maxQuantileAcc(ke);
         meanAcc(ke);
         }
    }
    allSpringsCalculated = (sumRecvBufferSize == 0);
  } while( !allSpringsCalculated );
  
  m_springTruncationData.minLimit = p_square_quantile(minQuantileAcc);
  m_springTruncationData.maxLimit = p_square_quantile(maxQuantileAcc);
  m_springTruncationData.mean     = mean(meanAcc);

  MPI_Bcast(&m_springTruncationData.minLimit, 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&m_springTruncationData.maxLimit, 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&m_springTruncationData.mean    , 1, MPI_DOUBLE, 0, communicator);

  CFLog(VERBOSE, "MeshFittingAlgorithm: Spring min limit: "<<m_springTruncationData.minLimit 
              << "; Spring max limit: "  <<  m_springTruncationData.maxLimit 
              << "; Spring mean Value: " << m_springTruncationData.mean << "\n");
}

////////////////////////////////////////////////////////////////////////////////////
CFreal MeshFittingAlgorithmFRQ2Quads::computeSpringConstantTruncation(const Framework::Node* const firstNode,
						     const Framework::Node* const secondNode) 
  {
  CFAUTOTRACE;
  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  //if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
    const CFreal firstNodeValue  = nodalStates[firstNode->getLocalID()] [m_monitorVarID];
    const CFreal secondNodeValue = nodalStates[secondNode->getLocalID()][m_monitorVarID];
  //}
  return std::abs(secondNodeValue - firstNodeValue);
}

//////////////////////////////////////////////////////////////////////////////
CFuint MeshFittingAlgorithmFRQ2Quads::getFaceInCommonID(const Framework::Node* const firstNode,
						     const Framework::Node* const secondNode){
  
  typedef CFMultiMap<CFuint, CFuint> MapNodeFace;
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  typedef MapNodeNode::MapIterator mapItNode;
   
  typedef MapNodeFace::MapIterator mapIt;
  typedef MapNodeFace::MapIterator mapItN;

  bool found = false;
  std::pair<mapIt,mapIt > ite=m_mapNodeFace1.find(firstNode->getLocalID(), found);
  cf_assert(found);

   
  bool foundN = false;
  std::pair<mapItN,mapItN > iteN=m_mapNodeFace1.find(secondNode->getLocalID(), foundN);
  cf_assert(foundN);

  CFuint faceInCommonID;

  for (mapIt itE = ite.first; itE != ite.second; ++itE) {
	      	   for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
			if (itE->second == itNe->second){
				faceInCommonID = itNe->second;
			}
		    }
  }
  
  return faceInCommonID;
  
}

///////////////////////////////////////////////////////////////////////////////
CFreal MeshFittingAlgorithmFRQ2Quads::computeSpringConstant(const Framework::Node* const firstNode,
						     const Framework::Node* const secondNode) 
{

  CFAUTOTRACE;

  CFuint faceInCommonID = getFaceInCommonID(firstNode, secondNode); //get the face linking the two nodes together

  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();


  CFreal firstNodeValue = 0;
  CFreal secondNodeValue = 0;

  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();

  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");


  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;
  geoDataFace.idx = faceInCommonID;
  m_face = m_faceBuilder->buildGE();

  m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
  m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

  m_states[LEFT ] = m_cells[LEFT ]->getStates();
  m_states[RIGHT] = m_cells[RIGHT]->getStates();

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = m_frData->getInnerFacesStartIdxs();
  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;
    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));
    }
  }
 }

 std::vector< Framework::Node*  >* m_NodesLeft = m_cells[LEFT ] ->getNodes();
 std::vector< Framework::Node*  >* m_NodesRight = m_cells[RIGHT] ->getNodes();

  
 // reset states in nodal points
 for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	m_firstLeftState[iVar] = 0;
 	m_secondLeftState[iVar] = 0;
 	m_firstRightState[iVar] = 0;
 	m_secondRightState[iVar] = 0;
 }
  
 // compute states in nodal points
  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
 {        
	if((*m_NodesLeft)[iNode]->getLocalID() == firstNode->getLocalID()){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	      		m_firstLeftState[iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[LEFT])[iSol])[iVar];
	                  }      
	          }
         }
	if((*m_NodesRight)[iNode]->getLocalID() == firstNode->getLocalID()){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      		m_firstRightState[iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[RIGHT])[iSol])[iVar];
	                  }      
	          }
         }
	if((*m_NodesLeft)[iNode]->getLocalID() == secondNode->getLocalID()){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      		m_secondLeftState[iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[LEFT])[iSol])[iVar];
	                  }      
	          }
         }
	if((*m_NodesRight)[iNode]->getLocalID() == secondNode->getLocalID()){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      		m_secondRightState[iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[RIGHT])[iSol])[iVar];
	                  }      
	          }
         }
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
      m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
      
  }
  
  RealVector FluxVector(m_dim);
  RealVector NodalVector(m_dim);

  // check if the flux points orientation
  for (CFuint iDim=0 ; iDim<m_dim; iDim++){
	FluxVector[iDim] = m_flxPntCoords[m_nbrFaceFlxPnts-1][iDim] - m_flxPntCoords[0][iDim];
	NodalVector[iDim] = (*secondNode)[iDim] - (*firstNode)[iDim];
  }
  
  bool aligned = (MathTools::MathFunctions::innerProd(FluxVector, NodalVector) > 0);
  
  std::vector< std::vector< Framework::State* > > m_actualStatesFlxPnt;
	
  m_actualStatesFlxPnt.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_actualStatesFlxPnt[LEFT].push_back(new State());
    m_actualStatesFlxPnt[RIGHT].push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_actualStatesFlxPnt[LEFT][iFlx]->setLocalID(iFlx);
    m_actualStatesFlxPnt[RIGHT][iFlx]->setLocalID(iFlx);
  }
 
  // compute states in flux points
  if (!aligned){
	for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt){ 
		*(m_actualStatesFlxPnt[LEFT][iFlxPnt]) = *(m_cellStatesFlxPnt[LEFT][m_nbrFaceFlxPnts - iFlxPnt -1]);
		*(m_actualStatesFlxPnt[RIGHT][iFlxPnt]) = *(m_cellStatesFlxPnt[RIGHT][m_nbrFaceFlxPnts - iFlxPnt -1]);
	}
  }
  else{
	for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt){ 
		*(m_actualStatesFlxPnt[LEFT][iFlxPnt]) = *(m_cellStatesFlxPnt[LEFT][iFlxPnt]);
		*(m_actualStatesFlxPnt[RIGHT][iFlxPnt]) = *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]);
	}
  }

  // compute the stiffness between the first node and the first flux point
  CFreal stiffFirstFlux = std::abs((*(m_actualStatesFlxPnt[LEFT][0]))[m_monitorVarID] - m_firstLeftState[m_monitorVarID]) + 
                           std::abs((*(m_actualStatesFlxPnt[RIGHT][0]))[m_monitorVarID] - m_firstRightState[m_monitorVarID]);

  // compute the stiffnesses between all the flux pointss
  CFreal stiffInside = 0;
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts-1; ++iFlxPnt){ 
	stiffInside += std::abs((*(m_actualStatesFlxPnt[LEFT][iFlxPnt]))[m_monitorVarID] - (*(m_actualStatesFlxPnt[LEFT][iFlxPnt+1]))[m_monitorVarID]) + 
                        std::abs((*(m_actualStatesFlxPnt[RIGHT][iFlxPnt]))[m_monitorVarID] - (*(m_actualStatesFlxPnt[RIGHT][iFlxPnt+1]))[m_monitorVarID]);
  }
  
  // compute the stiffness between the second node and the last flux point
  CFreal stiffFluxSecond = std::abs((*(m_actualStatesFlxPnt[LEFT][m_nbrFaceFlxPnts-1]))[m_monitorVarID] - m_secondLeftState[m_monitorVarID]) + 
                           std::abs((*(m_actualStatesFlxPnt[RIGHT][m_nbrFaceFlxPnts-1]))[m_monitorVarID] - m_secondRightState[m_monitorVarID]);

  CFreal stiffness = pow(stiffFirstFlux + stiffInside + stiffFluxSecond, 1);

  m_faceBuilder->releaseGE(); 

  return stiffness;
}

/////////////////////////////////////////////////////////////
CFreal MeshFittingAlgorithmFRQ2Quads::computeSpringConstantBoundary(const Framework::Node* const firstNode,
						     const Framework::Node* const secondNode) 
  {

  CFAUTOTRACE;

  CFuint faceInCommonID;
  CFuint commonTRS;
  
  typedef CFMultiMap<CFuint, CFuint> MapNodeFace;
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  typedef CFMultiMap<CFuint, CFuint> MapNodeTRS;

  typedef MapNodeNode::MapIterator mapItNode;
   
  typedef MapNodeFace::MapIterator mapIt;
  typedef MapNodeFace::MapIterator mapItN;

  typedef MapNodeTRS::MapIterator mapItTRS;
  typedef MapNodeTRS::MapIterator mapItNTRS;

  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  bool found = false;
  std::pair<mapIt,mapIt > ite=m_mapNodeFace1.find(firstNode->getLocalID(), found);
  cf_assert(found);

  CFreal firstNodeValue = 0;
  CFreal secondNodeValue = 0;


  bool foundN = false;
  std::pair<mapItN,mapItN > iteN=m_mapNodeFace1.find(secondNode->getLocalID(), foundN);
  cf_assert(foundN);

  for (mapIt itE = ite.first; itE != ite.second; ++itE) {
	      	   for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
			if (itE->second == itNe->second){
				faceInCommonID = itNe->second;
				//cout << faceInCommonID << endl;
			}
		    }
  }

  bool foundTRS = false;
  std::pair<mapItTRS,mapItTRS > itTRS =m_mapNodeTRS1.find(firstNode->getLocalID(), foundTRS);
  cf_assert(foundTRS);

  bool foundNTRS = false;
  std::pair<mapItNTRS,mapItNTRS > itNTRS =m_mapNodeTRS1.find(secondNode->getLocalID(), foundNTRS);
  cf_assert(foundNTRS);

  for (mapItTRS ittrs = itTRS.first; ittrs != itTRS.second; ++ittrs) {
	      	   for (mapItNTRS itntrs = itNTRS.first; itntrs != itNTRS.second; ++itntrs) {
			if (ittrs->second == itntrs->second){
				commonTRS = itntrs->second;
				//cout <<  "TRS" << commonTRS << endl;
			}
		    }
  }

  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();

  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  FaceToCellGEBuilder::GeoData& facesData = m_faceBuilder->getDataGE();
  m_faceBuilder->getDataGE().isBoundary = true;
  Common::SafePtr<Framework::TopologicalRegionSet> wallFaces =
  Framework::MeshDataStack::getActive()->getTrs( m_unlockedBoundaryTRSs[commonTRS] );
  facesData.facesTRS = wallFaces;
  map< std::string , vector< vector< CFuint > > >&
  bndFacesStartIdxsPerTRS = m_frData->getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[wallFaces->getName()];
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;
  // number of TRs
  const CFuint nbTRs = wallFaces->getNbTRs();
  //cout << "NBTRS" << nbTRs << endl;
  cf_assert(bndFacesStartIdxs.size() == nbTRs);
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = wallFaces;
  geoData.isBoundary = true;

  geoData.idx = faceInCommonID;
  m_face = m_faceBuilder->buildGE();

  std::vector<Framework::Node*>& faceNodes = *m_face->getNodes();

  // get the neighbouring cell
  m_intCell = m_face->getNeighborGeo(0);

  // get the states in the neighbouring cell
  m_cellStates = m_intCell->getStates();
  for (m_orient = 0; m_orient < nbOrients ; ++m_orient)
  {
   for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntsLocalCoords[iFlx] = (*m_allCellFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
      }
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
   {
    // reset the extrapolated states
    *(m_cellStatesFlxPntBnd[iFlxPnt]) = 0.0;
    
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    // extrapolate the states to current flx pnt
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

      *(m_cellStatesFlxPntBnd[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]));
      //cout << *(m_cellStatesFlxPntBnd[iFlxPnt]) << endl;
    }
  }
 }

 //cout << (*m_cellStates)[1]->getLocalID()  << endl;

 std::vector< Framework::Node*  >* m_Nodes = m_intCell->getNodes();
 //cout << "State" << ((*(m_cellStatesFlxPntBnd[0]))[0])<< endl;
 ///cout << "NODES size" << m_Nodes->size() << endl;

 std::vector< CFreal>  m_firstState;
 std::vector< CFreal>  m_secondState;

 m_firstState.resize(totalNbEqs);
 m_secondState.resize(totalNbEqs);

 // Resize in the setup only because very costly. Better to define variables in the hh and to reset to 0 when adding.
  
  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
 {        
	if((*m_Nodes)[iNode]->getLocalID() == firstNode->getLocalID()){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	      		m_firstState[iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
	                  }      
	          }
         }
	if((*m_Nodes)[iNode]->getLocalID() == secondNode->getLocalID()){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      		m_secondState[iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
	                  }      
	          }
         }
 }

 for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
      m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
  }
  
  RealVector FluxVector(m_dim);
  RealVector NodalVector(m_dim);

  for (CFuint iDim=0 ; iDim<m_dim; iDim++){
	FluxVector[iDim] = m_flxPntCoords[m_nbrFaceFlxPnts-1][iDim] - m_flxPntCoords[0][iDim];
	NodalVector[iDim] = (*secondNode)[iDim] - (*firstNode)[iDim];
  }
  
  bool aligned = (MathTools::MathFunctions::innerProd(FluxVector, NodalVector) > 0);
  std::vector< Framework::State* >m_actualStatesFlxPntBnd;


  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_actualStatesFlxPntBnd.push_back(new State());
    m_actualStatesFlxPntBnd.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_actualStatesFlxPntBnd[iFlx]->setLocalID(iFlx);
    m_actualStatesFlxPntBnd[iFlx]->setLocalID(iFlx);
  }
  
  if (!aligned){
	for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt){ 
		*(m_actualStatesFlxPntBnd[iFlxPnt]) = *(m_cellStatesFlxPntBnd[m_nbrFaceFlxPnts - iFlxPnt -1]);
	}
  }
  else{
	for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt){ 
		*(m_actualStatesFlxPntBnd[iFlxPnt]) = *(m_cellStatesFlxPntBnd[iFlxPnt]);
	}
  }

  
  CFreal stiffFirstFlux = std::abs((*(m_actualStatesFlxPntBnd[0]))[m_monitorVarID] - m_firstState[m_monitorVarID]);

  CFreal stiffInside = 0;
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts-1; ++iFlxPnt){ 
	stiffInside += std::abs((*(m_actualStatesFlxPntBnd[iFlxPnt]))[m_monitorVarID] - (*(m_actualStatesFlxPntBnd[iFlxPnt+1]))[m_monitorVarID]);
  }
  
  CFreal stiffFluxSecond = std::abs((*(m_actualStatesFlxPntBnd[m_nbrFaceFlxPnts-1]))[m_monitorVarID] - m_secondState[m_monitorVarID]);

  CFreal stiffness = 2*pow((stiffFirstFlux + stiffInside + stiffFluxSecond), 1);
  m_faceBuilder->releaseGE();

  firstNodeValue  = nodalStates[firstNode->getLocalID()] [m_monitorVarID];
  secondNodeValue = nodalStates[secondNode->getLocalID()][m_monitorVarID];

  return stiffness;
}
//////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithmFRQ2Quads::truncateSpringConstant(const CFreal springConstant){
  const CFreal maxLimit = m_springTruncationData.maxLimit;
  const CFreal minLimit = m_springTruncationData.minLimit;
  const CFreal mean     = m_springTruncationData.mean;

  const CFreal truncatedSpringConstant = std::max(std::min(maxLimit/mean, springConstant/mean), minLimit/mean );
  return truncatedSpringConstant;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::solveLinearSystem(){
  assembleLinearSystem();
  m_lss->solveSys();
	
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::resizeSystemSolverToNodalData(){
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint rhsSize = std::max(nodes.size()*totalNbEqs,states.size()*totalNbEqs);
  if (rhs.size()!=rhsSize) rhs.resize(rhsSize);
  rhs = 0.;
  jacobMatrix->resetToZeroEntries();
  //cout<<" -----------------state size " << states.size()<< endl;
  //cout<<" -----------------RHS size " << rhsSize<< endl;
  //cout<<" -----------------totalNbEqs size " << totalNbEqs<< endl;
  //cout<<" -----------------nodes size " << nodes.size()<< endl;

} 
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::resizeSystemSolverToStateData(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  if (rhs.size()/totalNbEqs != states.size()) rhs.resize(states.size()*totalNbEqs); 
}



//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::computeNodeStates()
{

    
  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  //Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  RealVector counter(nodalStates.size()); counter =0.;

  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    
    std::vector< Framework::Node*  >* m_cellNodes = currCell->getNodes();
    
    const CFuint nbNodes = m_cellNodes->size(); 


    // get the states in this cell
    std::vector< Framework::State* >* m_cellStates = currCell->getStates();



    for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
      {      
	// extrapolate the left and right states to the flx pnts
	for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	  {
	    for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      nodalStates[(*m_cellNodes)[iNode]->getLocalID()][iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
	      counter[(*m_cellNodes)[iNode]->getLocalID()]+=1.;
	    }      
	  }
      }
    m_cellBuilder->releaseGE();
    
  }
  for (CFuint i=0 ; i<nodalStates.size() ; ++i){
    nodalStates[i] = nodalStates[i]/counter[i];
    //cout<<" nodalStates   "  << nodalStates[i][0] <<"   " << nodalStates[i][1]  << "    " << nodalStates[i][2] <<"   " << nodalStates[i][3]<<endl;
    //cout<<" nodalStates   "  << nodalStates[i][0] <<"   " << nodalStates[i][1]  << "    " << nodalStates[i][2] <<"   " << nodalStates[i][3]<<endl;
  }
}
///////////////////////////////////////////////////////////////////

RealVector  MeshFittingAlgorithmFRQ2Quads::computeIntersection(const  Framework::Node* const  a,const  Framework::Node* const  c,
						      const  Framework::Node* const  b, const  Framework::Node* const  d){

  RealVector xyi(2); xyi = 0.;
  CFreal a1 = ((*c)[XX+1]- (*a)[XX+1])/((*c)[XX+0]-(*a)[XX+0]);
  CFreal a2 = ((*d)[XX+1]- (*b)[XX+1])/((*d)[XX+0]-(*b)[XX+0]);
  CFreal b1 =  (*a)[XX+1] - a1* (*a)[XX+0];
  CFreal b2 =  (*b)[XX+1] - a2* (*b)[XX+0];
  xyi[0] = (b1-b2)/(a2-a1);
  xyi[1] = a1* xyi[0]+b1;
  return xyi;
  
}


//////////////////////////////////////////////////////////////////////////////
  
CFreal MeshFittingAlgorithmFRQ2Quads::computeElementArea2dQuads(const  Framework::Node* const  a,const  Framework::Node* const  c,
						       const  Framework::Node* const  b, const  Framework::Node* const  d){

  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector vectorDC(3); vectorDC = 0.;
  RealVector vectorDA(3); vectorDA = 0.;
  RealVector vectorBC(3); vectorBC = 0.;
  RealVector vectorBA(3); vectorBA = 0.;
  RealVector vectorBD(3); vectorBD = 0.;
  RealVector vectorAC(3); vectorAC = 0.;

  RealVector res2(3);
  RealVector res1(3);
  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    vectorDC[iDim]= (*d)[XX+iDim]-(*c)[XX+iDim];
    vectorDA[iDim]= (*d)[XX+iDim]-(*a)[XX+iDim];
    vectorBC[iDim]= (*b)[XX+iDim]-(*c)[XX+iDim];
    vectorBA[iDim]= (*b)[XX+iDim]-(*a)[XX+iDim];
    vectorAC[iDim]= (*c)[XX+iDim]-(*a)[XX+iDim];
    vectorBD[iDim]= (*b)[XX+iDim]-(*d)[XX+iDim];

  }
  MathTools::MathFunctions::crossProd(vectorDC,vectorDA ,res1);
  MathTools::MathFunctions::crossProd(vectorBC,vectorBA ,res2);
  const CFreal  area =((res1.norm2()/2)+ (res2.norm2()/2));

  return area;
}

////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithmFRQ2Quads::computeConstantquads(const  RealVector xyi ,
						  const  Framework::Node* const  firstNode,
						  const  Framework::Node* const  secondNode,
						  const Framework::Node* const  thirdNode,
						  const Framework::Node* const  fourthNode,
						  CFreal elementArea){
  
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim(); 
  RealVector vector1(3); vector1 = 0.;
  RealVector vector2(3); vector2 = 0.;
  RealVector vector3(3); vector3 = 0.;
  RealVector vector4(3); vector4 = 0.;
  RealVector vector5(3); vector5 = 0.;
  RealVector vector6(3); vector6 = 0.;
  CFreal torsionConstant;
  RealVector res(3);
  RealVector res1(3);
  RealVector res2(3);
  //cout << "nbDims" << nbDims << endl;
  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    vector1[iDim]= (xyi)[iDim]-(*firstNode)[XX+iDim];
    vector2[iDim]= (xyi)[iDim]-(*secondNode)[XX+iDim];
    
    vector3[iDim]= (*thirdNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector4[iDim]= (*thirdNode)[XX+iDim]-(*secondNode)[XX+iDim];

    vector5[iDim]= (*fourthNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector6[iDim]= (*fourthNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }
  //cout << "oui" << endl;
  MathTools::MathFunctions::crossProd(vector1, vector2,res);
  const CFreal  area = (res.norm2()/2);
  //cout << "non" << endl;
  MathTools::MathFunctions::crossProd(vector3, vector4,res);
  const CFreal  area1 = (res.norm2()/2);

  MathTools::MathFunctions::crossProd(vector5, vector6,res);
  const CFreal  area2 = (res.norm2()/2);


  CFreal k = (area/elementArea)*(vector1.norm2()*vector2.norm2()*vector1.norm2()*vector2.norm2())/(4*area*area);

  CFreal k1 = (area1/elementArea)*(vector3.norm2()*vector3.norm2()*vector4.norm2()*vector4.norm2())/(4*area1*area1);

  CFreal k2 =(area2/elementArea)*(vector5.norm2()*vector5.norm2()*vector6.norm2()*vector6.norm2())/(4*area2*area2);
  // if( m_thetaMid ){
    // choose this option for a semi torsional spring analogy based on the middle angle only
  //torsionConstant = k;
  // }
  //else{
    // choose this option to stiffer the mesh or on a pave mesh
    torsionConstant =((std::max(std::max(k1,k2),k)));
    // }
  return torsionConstant;
}
      
//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::assembleLockedNode(const Framework::Node* node){
     Framework::DataHandle< CFreal > stiffness = socket_stiffness.getDataHandle();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;
  
  //cout << "xfirst " << (*node)[0] << "yfirst " << (*node)[1] << endl;
  //cout<< "LN: Global node ID   " << globalID << "   Local ID    " << node->getLocalID()<<endl;
   stiffness[node->getLocalID()]=100.;
  //cout << "before Jacob" << endl;
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, 1.);
  }
  //Right hand side
  //cout << "after Jacob" << endl;
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = (*node)[XX+iDim];
  }
  //cout << "after RHS" << endl;
}

//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::assembleinRegionNode2DQuads(const Framework::Node* node){
   //cout << "HERE1" << endl;
   Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  //cout << "NodeDist" << nodeDistance[node->getLocalID()] << endl;
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
  //cout << "NodeISAD" << nodeisAD[node->getLocalID()] << endl;
  Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  CFreal sumOffDiagonalValues = 0.;
  CFreal sum = 0.;
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  typedef MapNodeNode::MapIterator mapIt;
  typedef MapNodeNode::MapIterator mapItN;
  typedef MapNodeNode::MapIterator mapItS;
  bool found = false;
  std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
  cf_assert(found);
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
  //cout << "neigboring size" << neighboringNodes.size() << endl;
  
  std::vector<Framework::Node*>::const_iterator it;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
   	const  Framework::Node* neighborNode = *it;
	//cout << "xSEC" << (*neighborNode)[0] << "ySEC " << (*neighborNode)[1] << endl;
  }
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator itN1;
  std::vector<Framework::Node*>::const_iterator itSa;
  //cout << "HERE2" << endl;
  //cout << "xfirst " << (*node)[0] << "yfirst " << (*node)[1] << endl;
  //cout << "xsec " << (*secondNode)[0] << "ysec" << (*secondNode)[1] << endl;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){ // loop on the nodes 
    CFreal torsionConstant=0;
    std::vector<const Framework::Node*> sharedNodes;
    sharedNodes.clear();
    const  Framework::Node* neighborNode = *it;
    //cout << "HERE3" << endl;
    //cout << "xsecond " << (*neighborNode)[0] << "y second " << (*neighborNode)[1] << endl;
    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
	const Framework::Node* neighborNodeOfN = *itN;
         //cout << "xSECNNN" << (*neighborNodeOfN )[0] << "ySECNNN " << (*neighborNodeOfN )[1] << endl;
    }
    //cout << "neigboringN size" << neighboringNodesOfN.size() << endl;
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){ // loop on the neighbours of the neighbour node
      const Framework::Node* neighborNodeOfN = *itN;
      
      const std::vector<Framework::Node*>& neighboringNodes2 = m_edgeGraph.getNeighborNodesOfNode(node);
       for(itN1=neighboringNodes2.begin(); itN1 != neighboringNodes2.end(); ++itN1){
	const Framework::Node* neighborNodeNewloop = *itN1;
	if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
	  sharedNodes.push_back(neighborNodeNewloop);
           //cout << "x common " << (*neighborNodeNewloop)[0] << "y common " << (*neighborNodeNewloop)[1] << endl;
	  //cout << "sharedNodes size" << sharedNodes.size() << endl;
	}
      }
    }
    if (sharedNodes.size()==4){
      //cout << "HERE4" << endl;
      for (mapIt it = ite.first; it != ite.second; ++it) {
	for(CFuint i=0; i<sharedNodes.size() ; ++i){
	  if (sharedNodes[i]->getLocalID() == it->second){
	    bool foundN = false;
	    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
	    cf_assert(foundN);
	    bool foundS = false;
	    std::pair<mapItS,mapItS > iteS=m_mapNodeNode1.find(sharedNodes[i]->getLocalID(), foundS);
	    cf_assert(foundS);
	    for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) { // loop on the edge-connected neighbours nodes of the direct neighbour
	      for (mapItS itSe = iteS.first; itSe != iteS.second; ++itSe) {  // loop on the edge-connected neighbour nodes of the shared node
		if (itNe->second == itSe->second){
		  for(itSa=neighboringNodes.begin(); itSa != neighboringNodes.end(); ++itSa){
		    if((*itSa)->getLocalID() == itNe->second && (*itSa)->getLocalID()!= node->getLocalID()){ // neighbour of both the shared node and the neighbour 															but not of the original node
		      RealVector d_node_neighbor(nbDims); d_node_neighbor=0.;
		      RealVector d_neighbor_itSa(nbDims); d_neighbor_itSa=0.;
		      for (CFuint iDim=0; iDim<nbDims; ++iDim){
			d_node_neighbor[iDim] = (*node)[XX+iDim] - (*neighborNode)[XX+iDim];
			d_neighbor_itSa[iDim] = (*neighborNode)[XX+iDim] - (*(*itSa))[XX+iDim];
		      }
		      // Activate to put ST on the smallest edge and Linear spring on the longest edge
		      //if( d_node_neighbor.norm2() < d_neighbor_itSa.norm2() ){
		      RealVector xyi=computeIntersection(neighborNode, sharedNodes[i],node , *itSa);
		      CFreal elementArea=computeElementArea2dQuads(neighborNode, sharedNodes[i],node , *itSa);
		      torsionConstant+=computeConstantquads(xyi, node , neighborNode,  sharedNodes[i]   ,   *itSa  , elementArea);
		      //}
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      //cout << "torsionConstant" << torsionConstant << endl;
      const CFreal springConstant =computeSpringConstant(node,neighborNode);
      const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant); 
      //cout << "normalizedSpringConstant" << normalizedSpringConstant << endl;
      CFreal f = 1.;
      if (m_smoothSpringNetwork) {
      	//CFreal f = (7.-2.)/(0.02-0.01)*nodeDistance[node->getLocalID()]+ 1.2 - (7.-2.)/(0.02-0.01)*0.01;  //0.09
	// FB :test case dependent
      	 CFreal f = (5.-1.)/(0.0001-m_acceptableDistance)*nodeDistance[node->getLocalID()]+ 1. - (5.-1.)/(0.0001-m_acceptableDistance)*m_acceptableDistance;  // 5
      }
      const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
      const CFreal stiff=normalizedSpringConstant+(pow(torsionConstant,1.)*pow(normalizedSpringConstant,f));  
      //cout << "kstiff" << stiff << endl;
      //const CFreal stifNormalized = (stiff);
      //const CFreal stiff = normalizedSpringConstant;
      sumOffDiagonalValues +=(stiff) ;
      sum += 0.;
      stiffness[node->getLocalID()]=stiff;
      for(CFuint iDim=0; iDim<nbDims; ++iDim){
	jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,stiff);
      }
      
    }
  }
  
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue = -sumOffDiagonalValues-sum;
    const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
  }
  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue2 = sumOffDiagonalValues+sum;
    const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(diagValue2);
  }
 }

////////////////////////////////////////////////////////////////////////::
 void MeshFittingAlgorithmFRQ2Quads::assembleInnerNode(const Framework::Node* node){
   //SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
     //getTrs("InnerCells");
   //const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
   // Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
   //Framework::DataHandle<Framework::Node*, Framework::GLOBAL> neighborNodes = socket_nodes.getDataHandle();
     const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
     std::vector<Framework::Node*>::const_iterator itN;
 
   Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
   Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
   Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
   const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
   const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
   const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
   //2D quads
   if(true){
     CFreal sumOffDiagonalValues = 0.;
     CFreal sum = 0.;
     typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
     typedef MapNodeNode::MapIterator mapIt;
     typedef MapNodeNode::MapIterator mapItN;
     typedef MapNodeNode::MapIterator mapItS;
     bool found = false;
     std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
     cf_assert(found);    
      for(itN=neighboringNodes.begin(); itN != neighboringNodes.end(); ++itN){ 
       const  Framework::Node* neighborNode = *itN;
       for (mapIt it = ite.first; it != ite.second; ++it) {
            if (neighborNode->getLocalID() == it->second){
	      const CFreal springConstant =computeSpringConstant(node,neighborNode);
              CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
	     //normalizedSpringConstant = pow(normalizedSpringConstant, m_order);
	      //cout<<" normalized spring constant  "<< normalizedSpringConstant << endl;
              sumOffDiagonalValues +=normalizedSpringConstant ;
              const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
	      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
	      for(CFuint iDim=0; iDim<nbDims; ++iDim){
	         jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,normalizedSpringConstant);
              }

             //jacobMatrix->printToScreen();
	      //cout<< "couple node added  "<< rowGlobalID <<"  " << colGlobalID<< endl;
	      stiffness[node->getLocalID()]=normalizedSpringConstant;
	    }
	  }
       }
     for(CFuint iDim=0; iDim<nbDims; ++iDim){     
       const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;
       //cout<< "diagvalues  "<<globalID << endl;
       jacobMatrix->addValue(globalID+iDim, globalID+iDim, -sumOffDiagonalValues);
     }

     //Right hand side
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
       //cout<<" RHS Filling " << endl;
       rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues);
     }
     //jacobMatrix->printToScreen();

   }
 }

/////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::assembleMovingInBoundaryNode(const Framework::Node* node){
   
  bool blocked = false;
  Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  //Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  const RealVector& nodeNormal = (m_mapNodeIDNormal[node->getLocalID()]); // get the normal of each moving in boundary node
  
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator it;
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
  std::vector<const Framework::Node*> sharedNodes;
  for(itN=neighboringNodes.begin(); itN != neighboringNodes.end(); ++itN){
    Framework::Node* neighborNode = *itN;
    RealVector dist(nbDims);
    const bool neighborIsBoundary = m_boundaryNodes.find(neighborNode) != m_boundaryNodes.end();
    if (neighborIsBoundary){
    	for (CFuint iDim=0; iDim<nbDims; ++iDim){
		dist[iDim] = (*neighborNode)[XX+iDim]-(*node)[XX+iDim];
	}
	if (dist.norm2() < m_equilibriumSpringLength*1 && blocked == false){
		cout << "blocked" << endl;
		blocked = true;
	}
    }
  } 
  
  if (!blocked){
          //cout << "HERE" << endl;
  //Choose the dependent dimension and the free dimension 
  //by doing so we avoid problems with the pivoting of the LU decomposition
  
  //Dependent dim is the dimension with highest normal value
  	CFuint dependentDim = 0;
  	for (CFuint i=1; i<nbDims; ++i){
    		dependentDim = (std::abs(nodeNormal[i]) > std::abs(nodeNormal[dependentDim])) ? i : dependentDim ; // compare the different projected normals and determines the dependant one
  	}
  //Free dims are the others 
  	std::vector<CFuint> freeDims;
  	for (CFuint iDim=0; iDim<nbDims; ++iDim){
    		if (iDim!=dependentDim){
       			freeDims.push_back(iDim);
    		}
  	}

  //Dependent dimension, following the line equation a*y + b*x + c*z= d,
  // where y is the dep dim, x is the free dim, a and b are the normals and d is the y-intersect
  	const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
  	CFreal y_intersect = nodeNormal[dependentDim]*(*node)[XX+dependentDim]; // y*n_y
  	for(CFuint iFreeDim=0; iFreeDim<freeDims.size();++iFreeDim){
    		y_intersect += nodeNormal[freeDims[iFreeDim]]*(*node)[XX+freeDims[iFreeDim]]; // +x*n_x+z*n_z
  	}
  	jacobMatrix->addValue(globalID+dependentDim, globalID+dependentDim, nodeNormal[dependentDim]); // n_y
  	for(CFuint iFreeDim=0; iFreeDim<freeDims.size();++iFreeDim){
    		jacobMatrix->addValue(globalID+dependentDim, globalID+freeDims[iFreeDim], nodeNormal[freeDims[iFreeDim]]);// n_x and n_z 
  	}
  	rhs[node->getLocalID()*totalNbEqs+XX+dependentDim] = y_intersect; // x*n_x + y*n_y + z*n_z

  

//Free dimensions, acting like an inner Node with projected spring constants
	//cout << "Moving" << endl;
	CFreal sumOffDiagonalValues = 0.;
  	typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
     	typedef MapNodeNode::MapIterator mapIt;
     	bool found = false;
     	std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
     	cf_assert(found);  
	const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
         std::vector<Framework::Node*>::const_iterator itN;
  	for(itN=neighboringNodes.begin(); itN != neighboringNodes.end(); ++itN){
	   //cout << "neighbor" << endl;
	   Framework::Node* neighborNode = *itN;
	   for (mapIt it = ite.first; it != ite.second; ++it) {
               if (neighborNode->getLocalID() == it->second){
		const bool neighborIsBoundary = m_boundaryNodes.find(neighborNode) != m_boundaryNodes.end();
		CFreal springConstant = 0;
		if (neighborIsBoundary){
			//cout << " BLD" << endl;
			//cout << "xfirst " << (*node)[0] << "yfirst " << (*node)[1] << endl;
			//cout << "xsec " << (*neighborNode)[0] << "ysec" << (*neighborNode)[1] << endl;
			springConstant = computeSpringConstantBoundary(node, neighborNode); 
			//cout << "spring constant" << springConstant << endl;
		}
		else{
			//cout << "NO" << endl;
			//cout << "xfirst " << (*node)[0] << "yfirst " << (*node)[1] << endl;
			//cout << "xsec " << (*neighborNode)[0] << "ysec" << (*neighborNode)[1] << endl;
    			springConstant = computeSpringConstant(node, neighborNode); 
			//cout << "spring onstant" << springConstant << endl;
                  }// spring constant
		//cout << "after" << endl;
		//cout << "xfirst " << (*node)[0] << "yfirst " << (*node)[1] << endl;
        		//cout << "xsec " << (*neighborNode)[0] << "ysec" << (*neighborNode)[1] << endl;
    		RealVector springDirection(nbDims); // vector from the node to the neighbour
    		for (CFuint iDim=0; iDim<nbDims; ++iDim){
      			springDirection[iDim] = (*neighborNode)[XX+iDim] - (*node)[XX+iDim] ;
    		}
    		springDirection.normalize(); // normalized vector from the node to its neighbour
    		RealVector projectedSpringDirection = springDirection - 
      		MathTools::MathFunctions::innerProd(springDirection, nodeNormal)*nodeNormal;

    		const CFreal dotProduct = MathTools::MathFunctions::innerProd(projectedSpringDirection, springDirection);

    		const CFreal physicalSpringConstant =  truncateSpringConstant(springConstant)*std::abs(dotProduct);
		//const CFreal physicalSpringConstant =  truncateSpringConstant(springConstant);
    		CFreal normalizedSpringConstant =  physicalSpringConstant;
		//cout << "it is normalized" << endl;
		RealVector dist(nbDims);
                  for (CFuint iDim=0; iDim<nbDims; ++iDim){
		dist[iDim] = (*neighborNode)[XX+iDim]-(*node)[XX+iDim];
		}
		const CFreal stiff= normalizedSpringConstant;
          	stiffness[node->getLocalID()]=stiff;
		sumOffDiagonalValues += stiff;
		//cout << "it is stiff" << endl;
    		const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
    		const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
    		for(CFuint iFreeDim=0; iFreeDim<freeDims.size(); ++iFreeDim){
      			jacobMatrix->addValue(rowGlobalID+freeDims[iFreeDim], colGlobalID+freeDims[iFreeDim], stiff);
    		}
		//cout << "it is jacob" << endl;
	   }
          }
  	}
	//cout << "end of this node" << endl;
  	for(CFuint iFreeDim=0; iFreeDim<freeDims.size(); ++iFreeDim){
    		CFreal diagValue = -sumOffDiagonalValues;
    		jacobMatrix->addValue(globalID+freeDims[iFreeDim], globalID+freeDims[iFreeDim], diagValue);
    		rhs[node->getLocalID()*totalNbEqs+XX+freeDims[iFreeDim]] = m_equilibriumSpringLength*sumOffDiagonalValues;
  	}
   
}
if (blocked){
	assembleLockedNode(node);
   }	

}    
/////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::updateNodePositions () {
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::updateNodePositions()\n"); 
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
  Framework::DataHandle < CFreal > rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
   CFreal f=1.;
    if (nodes[iNode]->isParUpdatable()) {
      Framework::Node& currNode = *nodes[iNode];
      bool exit = false;
      if (m_smoothNodalDisp){
	// FB : test case dependent

	if( nodeDistance[nodes[iNode]->getLocalID()] < 0.00012  &&  insideRegion(nodes[iNode])==false ){ 				
	  f=(1./(0.00012 -m_acceptableDistance)) * nodeDistance[nodes[iNode]->getLocalID()] - (m_acceptableDistance/(0.00012  -m_acceptableDistance)); 
	}
      }
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
	currNode[XX+iDim] =  currNode[XX+iDim]*(1.-m_meshAcceleration*f) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration*f;
      }
      
    }
  }
  //synchronize Nodes
  nodes.beginSync();
  nodes.endSync();
}
/////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::assembleLinearSystem(){
  CFAUTOTRACE;
  //CFLog(VERBOSE, "MeshFittingAlgorithm::assembleLinearSystem()\n");
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
  Framework::DataHandle <CFreal> nodeDistance = socket_nodeDistance.getDataHandle();

  Common::SafePtr<Framework::TopologicalRegionSet> cells =
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
 //cout << "NB nodes " << nbNodes << endl;
 for (CFuint iNode = 0; iNode < nbNodes; ++iNode){
         //cout << "iNode" << iNode << endl;
	if (!nodes[iNode]->isParUpdatable()){ 
	  //do nothing
	}
  else{
	if(isNodeMovingInBoundary(nodes[iNode])){ //insideRegion(nodes[iNode])==false){
	  //cout << "moving" << endl;	
	  //cout << "xfirst " << (*nodes[iNode])[0] << "yfirst " << (*nodes[iNode])[1] << endl;
	  assembleMovingInBoundaryNode(nodes[iNode]);
	}	
	else{ 
	  if(isNodeLocked(nodes[iNode])){
	    //cout << "before locked" << endl;
	    assembleLockedNode(nodes[iNode]);
	    //cout << " after locked" << endl;	
		//cout << "xfirst " << (*nodes[iNode])[0] << "yfirst " << (*nodes[iNode])[1] << endl;
	  }
	  else{
	    if(insideRegion(nodes[iNode]) ){
	      if(nbElemTypes==3 && nbDims==2){
		//assembleinRegionNode2DTriag(nodes[iNode]);
	      }
	      if(nbElemTypes==4 && nbDims==2){
		//assembleLockedNode(nodes[iNode]);
                  //assembleinRegionNode2DQuads(nodes[iNode]);
		//cout << "before inner" << endl;
		assembleInnerNode(nodes[iNode]);
		//assembleInnerNode(nodes[iNode]);
		//cout << "after inner" << endl;
		//cout << "quads" << endl;	
		//cout << "xfirst " << (*nodes[iNode])[0] << "yfirst " << (*nodes[iNode])[1] << endl;
	      }
	      if(nbElemTypes==4 && nbDims==3){
		//assembleinRegionNode3DTet(nodes[iNode]);
	      }
	      if(nbElemTypes==8 && nbDims==3){
		//assembleinRegionNode3DHexa(nodes[iNode]);
	      }
	    }
	    else{
	      //cout << "inner" << endl;	
	      //cout << "xfirst " << (*nodes[iNode])[0] << "yfirst " << (*nodes[iNode])[1] << endl;
	      //cout << "before inner" << endl;
	      assembleInnerNode(nodes[iNode]);
	      //cout << "after inner" << endl;
	      //assembleinRegionNode3DHexa(nodes[iNode]);
	      //assembleinRegionNode2DQuads(nodes[iNode]);

	      //assembleinRegionNode2DQuads(nodes[iNode]);	
	      // assembleinRegionNode2DTriag(nodes[iNode]);
	    }	    
	 }
       }
     }
  }
}
/////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFRQ2Quads::isNodeLocked( Framework::Node* node)
{
  const bool isBoundary = m_boundaryNodes.find(node) != m_boundaryNodes.end();
  return ((isBoundary)) ;
}

//////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFRQ2Quads::insideRegion(Framework::Node* node)
{
  bool inRegion=false;
  Framework::Node& currNode = *node;
  const CFuint nodeID = currNode.getLocalID();
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
 // FB : uncomment isNodeLocked(node)==false f you are not using nodal movememt interpolation
  if  ((nodeisAD[nodeID] == true)){  //&& (isNodeLocked(node)==false)){
    inRegion=true;	    
  }
  return (inRegion);
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2Quads::findBoundaryNodes()
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::createConnectivity()" << "\n");
  
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  std::vector<CFuint> alreadyComputed;
  CFuint NBNODES = 0;
  for (CFuint iCell=0; iCell<nbCells; ++iCell) {
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    const std::vector<Framework::GeometricEntity*>* facesInCell = currCell->getNeighborGeos();
    Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryCell = 
      m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();


    const CFuint nbFaces = facesInCell->size();
    for (CFuint iFace=0; iFace<nbFaces; ++iFace) {
	//cout << "iFace" << iFace << endl;
      if((*m_isFaceOnBoundaryCell)[iFace]) {
	//cout << "iFace" << iFace << endl;
	std::vector<Framework::Node*>* faceNodes = (*facesInCell)[iFace]->getNodes();
	const CFuint nbFaceNodes = faceNodes->size();
	//cout << "numberNodalFaces" << nbFaceNodes << endl;
	bool done = false;
	for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {
	  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode((*faceNodes)[iNode]);
	  //cout << "neighbor size" << neighboringNodes.size() << endl;
	  //cout << "xFCELL " << (*(*faceNodes)[iNode])[0] << "yFCELL " << (*(*faceNodes)[iNode])[1] << endl;
	  for(CFuint jxx=0; jxx<alreadyComputed.size(); ++jxx){
	    if(alreadyComputed[jxx] ==  (*faceNodes)[iNode]->getLocalID()){
	      done = true;
	    }
	  }

	  if(done==false &&  neighboringNodes.size()<6 ){

	    	alreadyComputed.push_back((*faceNodes)[iNode]->getLocalID());
	    	m_boundaryNodes.insert((*faceNodes)[iNode]);
		//cout << "xFCELL " << (*(*faceNodes)[iNode])[0] << "yFCELL " << (*(*faceNodes)[iNode])[1] << endl;
		NBNODES +=1;
	  }
	}
      }
    }
    m_cellBuilder->releaseGE();
  }
  //cout << NBNODES << endl;
}


//////////////////////////////////////////////////////////////////////////////
   
void MeshFittingAlgorithmFRQ2Quads::computeMovingInBoundaryNodeNormals()
{

// FB : NOT Active YET
  CFAUTOTRACE;
  const CFuint nbDim = Framework::PhysicalModelStack::getActive()->getDim();
  //Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle(); 
  std::vector<std::set<CFuint> > mapNodeID2BoundaryFaceIDs = getMapMovingInBoundaryNodeID2FaceID();
  //Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle(); 
  std::set<Framework::Node*>::iterator it;
  RealVector normal1(nbDim), normal2(nbDim);
  //cout << "here" << endl;
  for(CFuint iNodeID=0; iNodeID<mapNodeID2BoundaryFaceIDs.size(); ++iNodeID){
    //Check if all normals of the connected faces are collinear up to a certain tolerance
    //that means check if the dot product of the normals are bellow that tolerance

    //cout << "here2" << endl;
    if(mapNodeID2BoundaryFaceIDs[iNodeID].size()>=nbDim){

      std::set<CFuint>& connectedFaceIds = mapNodeID2BoundaryFaceIDs[iNodeID];
      bool areNormalsCollinear  = true;
      std::set<CFuint>::iterator itFaceID1, itFaceID2;
      for(itFaceID1=connectedFaceIds.begin(); itFaceID1!=connectedFaceIds.end(); ++itFaceID1){
        for(itFaceID2=connectedFaceIds.begin(); itFaceID2!=connectedFaceIds.end(); ++itFaceID2){
          if(itFaceID1 != itFaceID2){
	   //cout << "faceDifferent" << endl;
            for(CFuint iDim=0; iDim<nbDim; ++iDim){
	      //cout<< "  (*itFaceID1)                "<<(*itFaceID1)<<endl;
	      //cout<< "  (*itFaceID1)*nbdims +idim   "<<(*itFaceID1)*nbDim+iDim<<endl;
	      //cout<< normals.size() << endl;
	      //cout << "start" << endl;
              normal1[iDim] = m_normalsAMR[(*itFaceID1)*nbDim+iDim];
              normal2[iDim] = m_normalsAMR[(*itFaceID2)*nbDim+iDim];
            }
            normal1.normalize(); normal2.normalize();
            const CFreal cosAngle = MathTools::MathFunctions::innerProd(normal1, normal2);
	    //cout<< " cosine angle    " << cosAngle<< endl;
            areNormalsCollinear &= (1.- std::abs(cosAngle)< 1e-5);
          }
        }
      }
      if(areNormalsCollinear){
        std::pair<CFuint, RealVector> addPair(0, RealVector(nbDim));
        addPair.first = iNodeID;
        //Compute Average Normal
        RealVector& averageNormal = addPair.second;
        std::set<CFuint>::iterator itFaceID;
        for(itFaceID=connectedFaceIds.begin(); itFaceID!=connectedFaceIds.end(); ++itFaceID){
          for(CFuint iDim=0; iDim<nbDim; ++iDim){
            const CFuint normalID = *itFaceID;
            averageNormal[iDim] += m_normalsAMR[normalID*nbDim+iDim];

          }
        }
        averageNormal.normalize();
        //CFLog(NOTICE,"averageNormal  "<<averageNormal<< "/n");
        m_mapNodeIDNormal.insert(addPair);
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
 
std::vector<std::set<CFuint> > MeshFittingAlgorithmFRQ2Quads::getMapMovingInBoundaryNodeID2FaceID()
{
  
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::computeNodes2BoundaryNormals()" << "\n");
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();
  std::vector<std::set<CFuint> > mapBoundaryNodeID2FaceIDs(nbNodes);
  //CFuint nbrNodes = 0;
  Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle(); 
  const CFuint nbDim = Framework::PhysicalModelStack::getActive()->getDim();
  //nodalStates[(*m_cellNodes)[iNode]->getLocalID()][iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
  for(CFuint iTrs=0; iTrs<m_unlockedBoundaryTRSs.size();++iTrs){
    FaceToCellGEBuilder::GeoData& facesData = m_faceBuilder->getDataGE();
    m_faceBuilder->getDataGE().isBoundary = true;
    Common::SafePtr<Framework::TopologicalRegionSet> wallFaces =
    Framework::MeshDataStack::getActive()->getTrs( m_unlockedBoundaryTRSs[iTrs] );
    facesData.facesTRS = wallFaces;
    const CFuint nbFaces = wallFaces->getLocalNbGeoEnts();
    //cout << "NBR BND Faces" << nbFaces << endl;
    for(CFuint i=0; i<nbFaces; ++i){
      facesData.idx = i;
      const Framework::GeometricEntity *const face = m_faceBuilder->buildGE();
      const std::vector<Framework::Node*> faceNodes = face->getNodes();
      const CFuint faceID = face->getID();
      const CFuint nbFaceNodes = faceNodes.size();
      for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {
	const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(faceNodes[iNode]);
	if(neighboringNodes.size()< 6){
          const CFuint nodeID = faceNodes[iNode]->getLocalID();
          //cout << "xfirst " << (*faceNodes[iNode])[0] << "yfirst " << (*faceNodes[iNode])[1] << endl;
	 //nbrNodes +=1;
          (mapBoundaryNodeID2FaceIDs[nodeID]).insert(faceID);
        }
      }
      //cout << "xfirst " << (*faceNodes[1])[0] << "yfirst " << (*faceNodes[1])[1] << endl;
      //cout << "xsec " << (*faceNodes[0])[0] << "ysec" << (*faceNodes[0])[1] << endl;
      CFreal deltaX = (*faceNodes[1])[0] - (*faceNodes[0])[0];
      CFreal deltaY = (*faceNodes[1])[1] - (*faceNodes[0])[1];
      //cout << "deltaX" << deltaX << endl;
      //cout << "deltaY" << deltaY << endl;
      m_normalsAMR[faceID*nbDim+0] = -deltaY;
      m_normalsAMR[faceID*nbDim+1] = deltaX;
      m_faceBuilder->releaseGE();
    }
  }
  //cout << "Nbr NODES" << nbrNodes << endl;
  return mapBoundaryNodeID2FaceIDs;

}

//////////////////////////////////////////////////////////////////////////////

bool MeshFittingAlgorithmFRQ2Quads::isNodeMovingInBoundary(Framework::Node* node){
  
  return (m_mapNodeIDNormal.find(node->getLocalID()) != m_mapNodeIDNormal.end());
}

/////////////////////////////////////////////// 
void MeshFittingAlgorithmFRQ2Quads::saveOldMeshProperties(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  cout << "begin" << endl;
  for(CFuint iState = 0; iState < states.size(); ++iState) {
    for (CFuint iVar = 0 ; iVar < nbEqs ; ++iVar){
      oldStates[nbEqs*states[iState]->getLocalID()+iVar]= (*states[iState])[iVar];
      //cout<< "saveOldMeshProperties::oldStates[nbEqs*states[iState]->getLocalID()+iVar]   "<<oldStates[nbEqs*states[iState]->getLocalID()+iVar]<<endl;
    }
    for (CFuint iDim=0; iDim< nbDims; ++iDim){ 
      RealVector& stateCoord =  states[iState]->getCoordinates();
      oldCoordinates[nbDims*states[iState]->getLocalID()+iDim] = stateCoord[iDim] ;
    }
  }
  cout << "end" << endl;
}

////////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2Quads::interpolateNewMeshProperties(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();;
  for(CFuint iState = 0; iState < states.size(); ++iState) {
    bool stopLooking =false;
    bool needInterpol =true;
    CFuint iStateOld = 0;
    CFreal invDistance = 0.;
    const RealVector& stateCoord =  states[iState]->getCoordinates();
    RealVector interpolatedStates(nbEqs); interpolatedStates=0.;
    while(iStateOld < states.size() && !stopLooking) {
      CFuint jxx = 0;
      for (CFuint iDim = 0; iDim < nbDims ; ++iDim){
	if(std::abs(oldCoordinates[nbDims*states[iStateOld]->getLocalID()+iDim]-stateCoord[iDim]) < 1e-8){
	  ++ jxx;
	}
      }
      if (jxx == nbDims){
	stopLooking = true;
	needInterpol = false;
	//cout<< " no Interpolation needed "<<states[iState]->getLocalID()<< endl;
	for(CFuint iVar = 0; iVar <nbEqs ; ++iVar) {
	  (*states[iState])[iVar] =  oldStates[nbEqs*states[iStateOld]->getLocalID()+iVar];
	}
      }
      ++ iStateOld;
    }

    if(needInterpol){
      CFuint closestState = 0;
      CFreal minDistance = MathTools::MathConsts::CFrealMax();
      for(CFuint iStateOld = 0; iStateOld < states.size(); ++iStateOld) {
	RealVector oldStateCoord(nbDims); oldStateCoord=0;
	for (CFuint iDim = 0; iDim < nbDims ; ++iDim){
	  oldStateCoord[iDim] = oldCoordinates[nbDims*states[iStateOld]->getLocalID()+iDim];
	}
	CFreal  newStateToOldState = MathTools::MathFunctions::getDistance(oldStateCoord , stateCoord );
	if (newStateToOldState < minDistance ){  
	  minDistance = newStateToOldState;
	  closestState = iStateOld;
	}
      }

     // cout<< " minDistance " << minDistance << endl;
      for(CFuint iVar = 0; iVar <nbEqs ; ++iVar) {
	//cout<< "  OLD(*states[iState])[iVar]  " <<  (*states[iState])[iVar]  << endl;
	//(*states[iState])[iVar] = interpolatedStates[iVar]/invDistance;
	(*states[iState])[iVar] =  oldStates[nbEqs*states[closestState]->getLocalID()+iVar];
	//cout<< "  NEW(*states[iState])[iVar]  " <<  (*states[iState])[iVar]  << endl;
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

 void MeshFittingAlgorithmFRQ2Quads::createGeneralConnectivityFR()
{ 
  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR(); 
   
  m_edgeGraphN.setNodeDataSocketFR(socket_nodes);
  m_edgeGraphN.computeConnectivityFR(); 
  }


//////////////////////////////////////////////////////////////////////////////
 
void MeshFittingAlgorithmFRQ2Quads::triggerRecomputeMeshData() {
  std::string msg;
  Common::SafePtr<Common::EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  const std::string ssname = Framework::SubSystemStatusStack::getCurrentName();   
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MESHADAPTER_AFTERMESHUPDATE"), msg );
  /*cout << "A0" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_UNSETUP"), msg );
  //cout << "A1" << endl;
 // event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_UNPLUGSOCKETS"), msg );
  //cout << "A2" << endl;
  // event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_DESTROYSUBSYSTEM"), msg ); 
 // cout << "A3" << endl;
  // event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_BUILDSUBSYSTEM"), msg );     
 //cout << "A4" << endl;
  // event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_PLUGSOCKETS"), msg );
  cout << "A5" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_BUILDMESHDATA"), msg );
  cout << "A6" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_SETUP"), msg );
  cout << "A7" << endl;*/
}
    }/// Namespace FR
  /// Namespace COOLFluiD
  
}
