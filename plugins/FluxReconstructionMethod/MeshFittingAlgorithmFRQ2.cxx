#include "Common/PE.hh"
#include "Common/EventHandler.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

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
#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolver.hh" 

#include "FluxReconstructionMethod/FluxReconstruction.hh"

#include "FluxReconstructionMethod/MeshFittingAlgorithmFRQ2.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "MeshTools/ComputeWallDistanceVectorFRMPI.hh"
#include "MeshTools/MeshToolsFR.hh"

#include "Common/CFMultiMap.hh"
#include <math.h>
#include <cmath>
#include <fstream>


///////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

///////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
        
///////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshFittingAlgorithmFRQ2, 
		      DataProcessingData,
		      FluxReconstructionModule>
MeshFittingAlgorithmQ2FRFluxReconstructionProvider("MeshFittingAlgorithmFRQ2");

///////////////////////////////////////////////////////////////////


void MeshFittingAlgorithmFRQ2::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("minPercentile","Percentile for minimum spring value");
  options.addConfigOption< CFreal >("maxPercentile","Percentile for maximum spring value");
  options.addConfigOption< CFreal >("meshAcceleration","How fast the mesh moves in mesh steps");
  options.addConfigOption< CFuint >("monitorVarID","Monitor variable ID (from State) for mesh adaptation");
  options.addConfigOption< CFuint >("monitorPhysVarID","Monitor physical variable ID (from physical data) for mesh adaptation");
  options.addConfigOption< CFreal >("equilibriumSpringLength","Length of spring for equilibrium");
  options.addConfigOption< CFreal >("ratioBoundaryToInnerEquilibriumSpringLength","ratio between the equilibrium length of a Boundary spring to an Inner spring");
  options.addConfigOption< std::vector<std::string> >("unlockedBoundaryTRSs","TRS's to be unlocked");
  options.addConfigOption< CFreal >("AcceptableDistanceQ2","Distance from user-defined boundary");
  options.addConfigOption< bool   >("ThetaMid","Semi torsional Sping analogy for 2D quadrilateral mesh based on the middle egde-facing  angle (true) or the 3 edge-facing angles (false) ");
  options.addConfigOption< bool   >("InterpolateState"," State Interplation to dissociate the nodal movement and the solution in each CC. ");
  options.addConfigOption< bool >("smoothSpringNetwork","smooth the spring network");
  options.addConfigOption< bool >("smoothNodalDisp","smooth the nodal displacememt");
}

///////////////////////////////////////////////////////////////////

MeshFittingAlgorithmFRQ2::MeshFittingAlgorithmFRQ2(const std::string& name) :

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
  m_secondRightState(),
  m_firstState(),
  m_secondState(),
  m_vecNodeCoords(),
  m_nodalCellStates(),
  m_NodeCoords()
  //m_faceTRSBuilder()
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

  m_acceptableDistanceQ2 = 0.;
  this->setParameter("AcceptableDistanceQ2",&m_acceptableDistanceQ2);

  m_thetaMid = true;
  this->setParameter("ThetaMid",&m_thetaMid);

  m_interpolateState = false;
  this->setParameter("InterpolateState",&m_interpolateState);

  m_smoothSpringNetwork = false;
  this->setParameter("smoothSpringNetwork",&m_smoothSpringNetwork);

  m_smoothNodalDisp =  true; // true;
  this->setParameter("smoothNodalDisp",&m_smoothNodalDisp);
}

//////////////////////////////////////////////////////////////////////////////

MeshFittingAlgorithmFRQ2::~MeshFittingAlgorithmFRQ2()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

 Framework::DataProcessingCom::configure(args);
  
  cf_assert(m_minPercentile >= 0.);
  cf_assert(m_minPercentile < m_maxPercentile);
  cf_assert(m_maxPercentile <= 1.);
  cf_assert(m_meshAcceleration > 0. && m_meshAcceleration < 1.);
  cf_assert(m_acceptableDistanceQ2 >= 0);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
MeshFittingAlgorithmFRQ2::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_normals);
  result.push_back(&socket_rhs);
  result.push_back(&socket_nstates);
  result.push_back(&socket_wallDistance);
  result.push_back(&socket_nodeisAD);  
  result.push_back(&socket_nodeDistance);
  return result;
}
//////////////////////////////////////////////////////////////////////////////
std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
MeshFittingAlgorithmFRQ2::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_stiffness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::setup()
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

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  //////////cout<<" nodalstate size()   "<<nodalStates.size()<<endl;
  //std::vector<CFuint > m_nbOfNeighborCellsToaNode;
  m_nbOfNeighborCellsToaNode.resize(socket_nodes.getDataHandle().size());

  nbOfConnectedFaces.resize(socket_nodes.getDataHandle().size());

    // AL: this might be useless ... (@see MeshRigidMove/StdSetup.cxx)
  SubSystemStatusStack::getActive()->setMovingMesh(true);
   
  //m_lss = getMethodData().getLinearSystemSolver()[0];
  const std::string name = getMethodData().getNamespace();
  m_lss = getMethodData().getCollaborator<LinearSystemSolver>(name);  

  CFLog(VERBOSE, "MeshFittingAlgorithmFRQ2::setup() -----=> LSS is " << m_lss->getName() << "\n");
  
  Common::SafePtr<Framework::SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();

  Common::SafePtr<FluxReconstructionSolver> frsolver = spaceMethod.d_castTo<FluxReconstructionSolver>(); //# Change here
  cf_assert(frsolver.isNotNull());
  m_frData = frsolver->getData();
  

  vector< FluxReconstructionElementData* >& frLocalData = m_frData->getFRLocalData();

 
  getMethodData().getUpdateVarSet()->setup();

  if (m_monitorVarID == std::numeric_limits<CFuint>::max() || 
    m_monitorVarID > PhysicalModelStack::getActive()->getNbEq()) {
    CFLog(WARN, "MeshFittingAlgorithmFRQ2::setup() => monitorVarID not specified or invalid: will be set to 0\n");
    m_monitorVarID = 0;
  }
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);

  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();

  m_order = static_cast<CFuint>(order);

  m_cellBuilder = m_frData->getCellBuilder();
  m_geoBuilder =  m_frData->getStdTrsGeoBuilder();
  m_faceBuilder = m_frData->getFaceBuilder();
  m_faceBuilder1 = m_frData->getFaceBuilder();



  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR();

  std::vector<CFuint >  nbOfNeighborCellsToaNode = m_edgeGraph.getNbOfNeighborCells();
  m_nbOfNeighborCellsToaNode = nbOfNeighborCellsToaNode;

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
  createGeneralConnectivityFR();    // New
  nbOfConnectedFacesToaNode();


  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();


  createNodeFaceConnectivity();



  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle(); 
  CFuint nbNodes = nodes.size();
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  
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
  
  //////////cout << "Total FACES" <<  nbTotalFaces << endl;
  m_normalsAMR.resize(nbTotalFaces*m_dim);
  m_orientsAMR.resize(nbTotalFaces*m_dim);
  findBoundaryNodes();

  for (CFuint root = 0; root < m_nbProc; ++root) {
  	computeMovingInBoundaryNodeNormals();  
  }
    createNodalConnectivity();

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

  m_firstState.resize(m_nbrEqs);
  m_secondState.resize(m_nbrEqs);


 m_vecNodeCoords.resize(6);

 m_vecNodeCoords[0].resize(2);
 m_vecNodeCoords[0][0] =0.;
 m_vecNodeCoords[0][1] =0.;

 m_vecNodeCoords[1].resize(2);
 m_vecNodeCoords[1][0] = 1.;
 m_vecNodeCoords[1][1] = 0.;

 m_vecNodeCoords[2].resize(2);
 m_vecNodeCoords[2][0] = 0.;
 m_vecNodeCoords[2][1] = 1.;

 m_vecNodeCoords[3].resize(2);
 m_vecNodeCoords[3][0] = .5;
 m_vecNodeCoords[3][1] = 0.;

 m_vecNodeCoords[4].resize(2);
 m_vecNodeCoords[4][0] = 0.5;
 m_vecNodeCoords[4][1] = 0.5;

 m_vecNodeCoords[5].resize(2);
 m_vecNodeCoords[5][0] = 0.;
 m_vecNodeCoords[5][1] = 0.5;


  /*
  m_vecNodeCoords[0].resize(2);
  m_vecNodeCoords[0][0] =-1.;
  m_vecNodeCoords[0][1] =-1.;

  m_vecNodeCoords[1].resize(2);
  m_vecNodeCoords[1][0] =-1.;
  m_vecNodeCoords[1][1] = 1.;
    
  m_vecNodeCoords[2].resize(2);
  m_vecNodeCoords[2][0] = 1.;
  m_vecNodeCoords[2][1] = 1.;

  m_vecNodeCoords[3].resize(2);
  m_vecNodeCoords[3][0] = 1.;
  m_vecNodeCoords[3][1] =-1.;

  m_vecNodeCoords[4].resize(2);
  m_vecNodeCoords[4][0] = 0.;
  m_vecNodeCoords[4][1] =-1.;

  m_vecNodeCoords[5].resize(2);
  m_vecNodeCoords[5][0] = 1.;
  m_vecNodeCoords[5][1] = 0.;

  m_vecNodeCoords[6].resize(2);
  m_vecNodeCoords[6][0] = 0.;
  m_vecNodeCoords[6][1] = 1.;

  m_vecNodeCoords[7].resize(2);
  m_vecNodeCoords[7][0] =-1.;
  m_vecNodeCoords[7][1] = 0.;

  m_vecNodeCoords[8].resize(2);
  m_vecNodeCoords[8][0] = 0.;
  m_vecNodeCoords[8][1] = 0.;
*/


  m_solPolyValsAtNodes = frLocalData[0]->getSolPolyValsAtNode(m_vecNodeCoords);
 //   m_solPolyValsAtNodes = frLocalData[0]->getCoefSolPolyInNodes();

  /*////////cout << m_vecNodeCoords[0] << endl;
  ////////cout << m_vecNodeCoords[1] << endl;
  ////////cout << m_vecNodeCoords[2] << endl;
  ////////cout << m_vecNodeCoords[3] << endl;
  ////////cout << m_vecNodeCoords[4] << endl;
  ////////cout << m_vecNodeCoords[5] << endl;
  ////////cout << m_vecNodeCoords[6] << endl;
  ////////cout << m_vecNodeCoords[7] << endl;
  ////////cout << m_vecNodeCoords[8] << endl;*/

  ////////cout << m_vecNodeCoords.size() << endl;

  m_nbrNodesElem = m_solPolyValsAtNodes.size();

  std::vector< std::vector< Framework::State* > > m_edgeStiffness;

  DistanceFirstNodeFlux.resize(m_nbrFaceFlxPnts);

}

//////////////////////////////////////////////////////////////////////////////


void MeshFittingAlgorithmFRQ2::unsetup()
{
  CFAUTOTRACE;
  Framework::DataProcessingCom::unsetup();
}
  

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::createNodeFaceConnectivity()
{ 
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
  Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();

  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle(); 
  CFuint nbNodes = nodes.size();
  

  // Connectivity information 2D quadrilateral
  CFuint nbPairsNodeNode = 20000; 

  typedef CFMultiMap<CFuint, CFuint > MapFaceNode;
  Common::CFMultiMap<CFuint,CFuint>  m_mapFaceNode(20000);
  typedef MapFaceNode::MapIterator mapFaceIt;

  typedef CFMultiMap<CFuint, CFuint > MapNodeFace;
  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeFace(20000);
  typedef MapNodeFace::MapIterator mapNodeIt;

  typedef CFMultiMap<CFuint, CFuint > MapNodeTRS;
  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeTRS(20000);
  typedef MapNodeTRS::MapIterator mapNodeTRSIt;

  typedef CFMultiMap<CFuint, CFuint > MapNodeFaceTRS;
  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeFaceTRS(200000);
  typedef MapNodeFace::MapIterator mapNodeFaceTRSIt;  

  typedef CFMultiMap<CFuint, CFuint > MapNodeCell;
  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeCell(20000);
  typedef MapNodeCell::MapIterator mapNodeCellIt;
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  m_nodalCellStates.resize(nbCells);

  for (CFuint iN = 0; iN < nbNodes; ++iN)
  {
    for (CFuint i = 0; i<nbCells; i++){
    	m_nodalCellStates[i].push_back(new State());
    }
  }

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
        //cout<<"face ID inner "<< faceID<<endl;
      }
      m_mapFaceNode.sortKeys();
      m_mapNodeFace.sortKeys();
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
			      m_mapNodeFaceTRS.insert(faceNodes[iNodeC]->getLocalID(), faceID);
			      m_mapNodeTRS.insert(faceNodes[iNodeC]->getLocalID(),iTrs);
	        }
	        m_mapFaceNode.sortKeys();
	        m_mapNodeFace.sortKeys();
          m_mapNodeFaceTRS.sortKeys();
	        m_mapNodeTRS.sortKeys();
	        m_faceBuilder->releaseGE();
	      }
	  }
     }
  }
  m_mapNodeFace1 = m_mapNodeFace;
  m_mapFaceNode1 = m_mapFaceNode;
  m_mapNodeTRS1 = m_mapNodeTRS;
  m_mapNodeFaceTRS1= m_mapNodeFaceTRS; 
	      CFLog(INFO,"m_mapNodeFaceTRS1: " << m_mapNodeTRS1.size() << "\n");

  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR();
}
/////////////////////////////////////////////////////////////////////


 void MeshFittingAlgorithmFRQ2::createNodalConnectivity()
{ 
  /////////////////////////
  //2D quadrilateral mesh//
  /////////////////////////
            //////////cout<<" Start Connectvity "<<endl;
CFuint test_var =0;

  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle(); 

  // Connectivity information 2D quadrilateral
  CFuint nbPairsNodeNode = 20000; 
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  Common::CFMultiMap<CFuint, CFuint>  m_mapNodeNode(nbPairsNodeNode);
  
  std::multimap<CFuint, CFuint>  mapNodeNode;
  std::vector< bool > nodeDone; nodeDone.resize(nodes.size());
  std::vector< CFuint > BCnodeDone; BCnodeDone.resize(nodes.size());  // count conectivity to each node 


  for(CFuint i=0; i<BCnodeDone.size(); ++i){
    BCnodeDone[i] = 0;
  }

    for(CFuint i=0; i<nodeDone.size(); ++i){
    nodeDone[i] = false;
  }



  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  bool coupleDone [nodes.size()][nodes.size()];
  for (CFuint k= 0 ; k<nodes.size() ; ++k){
    for (CFuint j= 0 ; j<nodes.size() ; ++j){
      coupleDone[k][j]=false;
    }
  }


  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    std::vector< Framework::Node*  >* m_cellNodes = currCell->getNodes();
    const CFuint nbNodes = m_cellNodes->size(); 
    const std::vector<Framework::GeometricEntity*  >& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size(); 
    for (CFuint iFace=0; iFace<nbFaces; ++iFace){
      std::vector<Framework::Node* >& faceNodes = *facesInCell[iFace]->getNodes();
      const CFuint nbNodesinF = faceNodes.size();
     
      for (CFuint iNode=0; iNode<nbNodesinF; ++iNode){
        Framework::Node& currNode = *(faceNodes)[iNode];
        if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] == 1 && isBoundaryNode(&currNode) && nbOfConnectedFaces[currNode.getLocalID()] == 1){ // Our node: Boundary Middle Point
          bool alreadyComputed = false;
          for (CFuint iNode=0; iNode<nbNodesinF; ++iNode){
            if((faceNodes)[iNode]->getLocalID() != currNode.getLocalID()){
              m_mapNodeNode.insert(currNode.getLocalID(), (faceNodes)[iNode]->getLocalID()); // Adding the 2 end points 
              BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;

            }

            /*for(CFuint iNode=0; iNode< nbNodes; ++ iNode){
              if (isInSideCell((*m_cellNodes)[iNode]) && alreadyComputed == false){
                m_mapNodeNode.insert(currNode.getLocalID(),(*m_cellNodes)[iNode]->getLocalID()); // Adding the Center point 
                BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
                alreadyComputed = true;
              }
            }*/
          }
          nodeDone[currNode.getLocalID()] = true;
        }




        if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] >= 1 && isBoundaryNode(&currNode) && nbOfConnectedFaces[currNode.getLocalID()] != 1){ // Our node: Boundary intersection Point
          RealVector distance(nbDims); distance=MathTools::MathConsts::CFrealMax();
          CFreal minDistance  = MathTools::MathConsts::CFrealMax();
          CFuint closestNode;
          for (CFuint iNode=0; iNode<nbNodesinF; ++iNode){
            if((faceNodes)[iNode]->getLocalID() != currNode.getLocalID()){
              Framework::Node& n_currNode = *(faceNodes)[iNode];
              for (CFuint iDim = 0; iDim <nbDims;  ++iDim){
                distance[iDim] = (currNode)[XX+iDim] - (n_currNode)[XX+iDim];
              }
              if (distance.norm2() < minDistance){
                minDistance = distance.norm2();
                closestNode = (faceNodes)[iNode]->getLocalID();
              }
            }
          }
          if(coupleDone[currNode.getLocalID()][closestNode]==false){
            m_mapNodeNode.insert(currNode.getLocalID(),closestNode); // Adding the middle point either on left or right or down/up
            nodeDone[currNode.getLocalID()] = true;
            coupleDone[currNode.getLocalID()][closestNode]=true;
                         //cout<< " HERE BCNode 2 ------------------------"<<endl;;

            BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
          }
        }

     /* if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] == 1 && isBoundaryNode(&currNode) && nbOfConnectedFaces[currNode.getLocalID()] == 2){ // Our node: Corner Point
          RealVector distance(nbDims); distance=MathTools::MathConsts::CFrealMax();
          CFreal minDistance  = MathTools::MathConsts::CFrealMax();
          CFuint closestNode;
          for (CFuint iNode=0; iNode<nbNodesinF; ++iNode){
            if((faceNodes)[iNode]->getLocalID() != currNode.getLocalID()){
              Framework::Node& n_currNode = *(faceNodes)[iNode];
              for (CFuint iDim = 0; iDim <nbDims;  ++iDim){
                distance[iDim] = (currNode)[XX+iDim] - (n_currNode)[XX+iDim];
              }
              if (distance.norm2() < minDistance){
                minDistance = distance.norm2();
                closestNode = (faceNodes)[iNode]->getLocalID();
              }
            }
          }
          m_mapNodeNode.insert(currNode.getLocalID(),closestNode); // Adding the middle point either on left or right or down/up
          nodeDone[currNode.getLocalID()] = true;
          BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
        }*/

        if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] == 2 && isBoundaryNode(&currNode)==false   ){ // Our node: Middle point between 2 cells    && nodeDone[currNode.getLocalID()] == false        
          //&& nodeDone[currNode.getLocalID()] == false
        /*    CFuint count=0;

          const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(&currNode);
          std::vector<Framework::Node*>::const_iterator itN;
          //////cout<<"ICIC ----------------"<<neighboringNodes.size()<<endl;
          for(itN=neighboringNodes.begin(); itN != neighboringNodes.end(); ++itN){ 
            const Framework::Node* neighborNode = *itN;
            m_mapNodeNode.insert(currNode.getLocalID(), neighborNode->getLocalID());
            count+=1;
          }*/
           for (CFuint iNode=0; iNode<nbNodesinF; ++iNode){
             if((faceNodes)[iNode]->getLocalID() != currNode.getLocalID()   && coupleDone[currNode.getLocalID()][(faceNodes)[iNode]->getLocalID()]==false){
               m_mapNodeNode.insert(currNode.getLocalID(), (faceNodes)[iNode]->getLocalID()); // Adding the 2 end points 
               coupleDone[currNode.getLocalID()][(faceNodes)[iNode]->getLocalID()]=true;  // Cheking if they are done
             }
           }
            
        ////cout<<BCnodeDone[currNode.getLocalID()]<<endl;

          nodeDone[currNode.getLocalID()] = true;
        }




        if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] >= 4 && isBoundaryNode(&currNode)==false){ // Our node: intersection between 2 edges
          
          RealVector distance(nbDims); distance=MathTools::MathConsts::CFrealMax();
          CFreal minDistance  = MathTools::MathConsts::CFrealMax();
          CFuint closestNode;
          for (CFuint iNodeF=0; iNodeF<nbNodesinF; ++iNodeF){
            if((faceNodes)[iNodeF]->getLocalID() != currNode.getLocalID()){
              Framework::Node& n_currNode = *(faceNodes)[iNodeF];
              for (CFuint iDim = 0; iDim <nbDims;  ++iDim){
                distance[iDim] = (currNode)[XX+iDim] - (n_currNode)[XX+iDim];
              }
              if (distance.norm2() < minDistance){
                minDistance = distance.norm2();
                closestNode = (faceNodes)[iNodeF]->getLocalID();
              }
            }
          }
          if(coupleDone[currNode.getLocalID()][closestNode]==false){
          m_mapNodeNode.insert(currNode.getLocalID(),closestNode); // Adding the middle point either on left or right or down/up
          nodeDone[currNode.getLocalID()] = true;
          coupleDone[currNode.getLocalID()][closestNode]=true;

          }
        }
      }
    }
    m_cellBuilder->releaseGE();
  }
  /*for (CFuint iNode=0; iNode<nodes.size(); ++iNode){
    Framework::Node& currNode = *nodes[iNode];
      if(nodeDone[currNode.getLocalID()]==false && isBoundaryNode(nodes[iNode])==false){ // Our node is a center node
        
        for (CFuint iCell=0; iCell<nbCells; ++iCell){
          geoData.idx = iCell;
          bool istheCellwiththeCurrentNode = false;
          Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
          std::vector< Framework::Node*  >& m_cellNodes = *currCell->getNodes();

          for(CFuint iNodeinCell=0; iNodeinCell<m_cellNodes.size(); ++iNodeinCell){
            if(m_cellNodes[iNodeinCell]->getLocalID() == currNode.getLocalID() &&  isBoundaryNode(m_cellNodes[iNodeinCell])==false ){
              istheCellwiththeCurrentNode = true;
            }
          }
          if(istheCellwiththeCurrentNode){

          const std::vector<Framework::GeometricEntity*  >& facesInCell = *currCell->getNeighborGeos();
          const CFuint nbFaces = facesInCell.size(); 
          for (CFuint iFace=0; iFace<nbFaces; ++iFace){
	          std::vector<Framework::Node* >& faceNodes = *facesInCell[iFace]->getNodes();
            const CFuint nbNodesinF = faceNodes.size();
            RealVector distance(nbDims); distance=MathTools::MathConsts::CFrealMax();
            CFreal minDistance  = MathTools::MathConsts::CFrealMax();
            CFuint closestNode;
            for (CFuint iNode=0; iNode<nbNodesinF; ++iNode){
              Framework::Node& neighborOfcurrNode = *(faceNodes)[iNode];  // nodes of the faces of our cell 
              for (CFuint iDim = 0; iDim <nbDims;  ++iDim){
                distance[iDim] = (currNode)[XX+iDim] - (neighborOfcurrNode)[XX+iDim];
              }
              if (distance.norm2() < minDistance){
                minDistance = distance.norm2();
                closestNode = (faceNodes)[iNode]->getLocalID();
              }
            }
          m_mapNodeNode.insert(currNode.getLocalID(),closestNode); // Adding the closest point within an edge to the center points
          BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
          }
          nodeDone[currNode.getLocalID()] = true;
          }
          m_cellBuilder->releaseGE();
          
        }
      }
  }*/

  m_mapNodeNode.sortKeys();
  //   }
  m_mapNodeNode1=m_mapNodeNode;
   ////cout<<"NEW PRINT HERE:---------------------------------------------------------------- "<<m_mapNodeNode1.size()<<endl;
	      CFLog(INFO,"Finish Connectvity ----------------" << "\n");
for (CFuint i=0; i<BCnodeDone.size(); i++) {
  //if (BCnodeDone[i]!=0){
    //cout<< i <<"----------->  "<< BCnodeDone[i]<<endl;
  //}
}
          //////////cout<<" Finish Connectvity "<<endl;
          ////////cout<<" test var "<< test_var<<endl;
}

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::execute()  /// needs to be changed 
{
  CFAUTOTRACE;
  CFLog(INFO, "MeshFittingAlgorithmFRQ2::execute() => start \n");

  determineIsNodeAD();

  resizeSystemSolverToNodalData();
    //CFLog(INFO, "here 1 \n");

  computeNodeStates();
    //CFLog(INFO, "here 2 \n");

  computeSpringTruncationData();
    //CFLog(INFO, "hi 2 \n");

  solveLinearSystem();
  // mesh node repositionning
    //CFLog(INFO, "hi 1 \n");

  updateNodePositions();

  resizeSystemSolverToStateData();

  triggerRecomputeMeshData();
  CFLog(INFO, "MeshFittingAlgorithmQ2::execute() => end \n");
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::determineIsNodeAD(){
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<bool> nodeisAD = socket_nodeisAD.getDataHandle();
  
  
   for (CFuint iNode = 0; iNode < nodes.size(); ++iNode){
	nodeisAD[nodes[iNode]->getLocalID()] = false;
	//////////cout << "Node Dist" << nodeDistance[nodes[iNode]->getLocalID()] << endl;
	if (nodeDistance[nodes[iNode]->getLocalID()] < m_acceptableDistanceQ2){
		nodeisAD[nodes[iNode]->getLocalID()] = false; // FB: change here after 
   	}
	//////////cout << "NODE IS AD" << nodeisAD[nodes[iNode]->getLocalID()] << endl;
  
   }
}

//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::computeSpringTruncationData() 
{
  //This method computes the average, min and max spring constants through a user defined quantile
  //The P^2 algorithm was used to compute the quantiles using the BOOST statistical accumulators
  //This method has very low memory requirements compared to a naive implementation but lower accuracy. 
  //The sum of accumulators are not defined so a parallel reduce can't be used
  //Instead, all the spring constants need to be gathered on process 0 to be accumulated there
  //In particular, the P^2 algorithm will be dependent on the ordering of the springs and thus,
  // on the number of processors and size of the send buffers

  CFAUTOTRACE;
  CFLogDebugMin( "MeshFittingAlgorithmQ2::computeSprings()" << "\n");
  
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle(); 

  using namespace boost::accumulators;
  typedef accumulator_set<CFreal, stats<tag::p_square_quantile> > accumulator_t;
  accumulator_t minQuantileAcc(quantile_probability = m_minPercentile);
  accumulator_t maxQuantileAcc(quantile_probability = m_maxPercentile);
  accumulator_set<CFreal, stats<tag::mean> > meanAcc;

  const std::string nsp = this->getMethodData().getNamespace();
  const int nbProcesses = Common::PE::GetPE().GetProcessorCount(nsp);
  const int processRank = Common::PE::GetPE().GetRank(nsp);
  MPI_Comm communicator = Common::PE::GetPE().GetCommunicator(nsp); 
  
  //Limit the size of the receiving buffer on processor 0
  const CFuint sizeSendBuffer = std::ceil(1000./static_cast<CFreal>(nbProcesses)) + 1;
  const CFuint sizeRecvBuffer = sizeSendBuffer*nbProcesses;
    //CFLog(INFO, "here 4 \n");
				      CFreal springConstant = 0;

  std::vector<CFreal> sendBuffer, recvBuffer; 
  sendBuffer.reserve(sizeSendBuffer);
  if(processRank == 0) recvBuffer.reserve(sizeRecvBuffer);
  
  SimpleEdgeGraphFRQ2::iterator it = m_edgeGraph.begin();
  bool allSpringsCalculated = false;
  do{
    //compute store the spring constants in a temporary buffer
    for(CFuint iBuffer=0; iBuffer<sizeSendBuffer && it != m_edgeGraph.end(); ++it){
      Framework::Node* node1 = (*it).firstNode;
      Framework::Node* node2 = (*it).secondNode;
      if(node1->isParUpdatable()){
	      typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
	      typedef MapNodeNode::MapIterator mapIt;
	      bool foundKPS = false;
            //CFLog(INFO, "here 5 \n");

	      std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node1->getLocalID(), foundKPS);
	      cf_assert(foundKPS); 
	      const bool nodeIsBoundary = (m_mapNodeIDNormal.find(node1->getLocalID()) != m_mapNodeIDNormal.end());
	      if(!isNodeLocked(node1) && !nodeIsBoundary){   
	     	  for (mapIt it = ite.first; it != ite.second; ++it) {
		   	    if (node2->getLocalID() == it->second){
				      //CFreal springConstant = 0;
    //CFLog(INFO, "here 6 \n");

				      springConstant = computeSpringConstant(node1,node2);
    //CFLog(INFO, "here 8 \n");
/*
	              if(springConstant > 1.e-14){
		 			        sendBuffer.push_back(springConstant);
		 			        ++iBuffer;
                       //CFLog(INFO, "here 10 \n");

                }*/
            }
          }
        }
	      if(isNodeLocked(node1) && nodeIsBoundary){
              //CFLog(INFO, "here 11 \n");

		      for (mapIt it = ite.first; it != ite.second; ++it) {
		   	    if (node2->getLocalID() == it->second){
				      const bool neighborIsBoundary = m_boundaryNodes.find(node2) != m_boundaryNodes.end();
                  //CFLog(INFO, "here 15 \n");

				      CFreal springConstant = 0;
				      if (neighborIsBoundary){
                CFuint commonTRS;
                bool sameTRS = false;

                typedef CFMultiMap<CFuint, CFuint> MapNodeTRS1;
                typedef MapNodeTRS1::MapIterator mapItTRS;
                typedef MapNodeTRS1::MapIterator mapItNTRS;

                bool foundTRS = false;
                std::pair<mapItTRS,mapItTRS > itTRS =m_mapNodeTRS1.find(node1->getLocalID(), foundTRS);
                cf_assert(foundTRS);

                bool foundNTRS = false;
                //CFLog(INFO, "node1: "<<(*node1)[XX] << "  "  <<  (*node1)[YY]<< "\n");

                //CFLog(INFO, "node2: "<<(*node2)[XX] << "  "  <<  (*node2)[YY]<< "\n");
                std::pair<mapItNTRS,mapItNTRS > itNTRS =m_mapNodeTRS1.find(node2->getLocalID(), foundNTRS);
                cf_assert(foundNTRS);
                for (mapItTRS ittrs = itTRS.first; ittrs != itTRS.second; ++ittrs) {
	                for (mapItNTRS itntrs = itNTRS.first; itntrs != itNTRS.second; ++itntrs) {
			              if (ittrs->second == itntrs->second){
				              commonTRS = itntrs->second;
                      sameTRS = true;
			              }
		              }
                }
                if (sameTRS){
					        springConstant =computeSpringConstant(node1, node2);
                }
                else{
                  springConstant =computeSpringConstant(node1, node2);
                }  
             
				      }
              
				      else {
					      springConstant =computeSpringConstant(node1, node2);
				      }
			      }
		     }
	      }
              if(springConstant > 1.e-4){
		 			        sendBuffer.push_back(springConstant);
		 			        ++iBuffer;
                       //CFLog(INFO, "here 10 \n");

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
    	////cout << "here 1--------" << endl;

    for(CFuint i=0; i< nbProcesses; ++i) {
      recvBufferDisp[i] = sumRecvBufferSize;
      sumRecvBufferSize += recvBufferCounts[i];
    }
    if(processRank == 0) recvBuffer.resize(sumRecvBufferSize); 
    //gather the temporary buffers in processor 0
    //to perform the statistical computations
        	////cout << "here 2--------" << endl;

    MPI_Gatherv(&sendBuffer[0], sendBufferSize, MPI_DOUBLE, &recvBuffer[0], 
        &recvBufferCounts[0], &recvBufferDisp[0], MPI_DOUBLE, 0, communicator); 
    sendBuffer.clear();  
    if(processRank == 0) {
	    for(CFuint i=0; i<recvBuffer.size(); ++i){
		    const CFreal ke = recvBuffer[i];
        minQuantileAcc(ke);
        maxQuantileAcc(ke);
        meanAcc(ke);
            	////cout << "here 3--------" << endl;

      }
    }
    allSpringsCalculated = (sumRecvBufferSize == 0);
  } while( !allSpringsCalculated );
  
  m_springTruncationData.minLimit = p_square_quantile(minQuantileAcc);
  m_springTruncationData.maxLimit = p_square_quantile(maxQuantileAcc);
  m_springTruncationData.mean     = mean(meanAcc);
    	////cout << "here 3--------" << endl;

  MPI_Bcast(&m_springTruncationData.minLimit, 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&m_springTruncationData.maxLimit, 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&m_springTruncationData.mean    , 1, MPI_DOUBLE, 0, communicator);

  CFLog(INFO, "MeshFittingAlgorithmQ2: Spring min limit: "<<m_springTruncationData.minLimit 
              << "; Spring max limit: "  <<  m_springTruncationData.maxLimit 
              << "; Spring mean Value: " << m_springTruncationData.mean << "\n");

              ////cout<<"Finished--------------------------------------------"<<endl;

}

/////////////////////////////////////////////////////////////////////////////////////////
CFuint MeshFittingAlgorithmFRQ2::getFaceInCommonID(const Framework::Node* const firstNode, 
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

//////////////////////////////////////////////////////////////////////////////

  CFreal MeshFittingAlgorithmFRQ2::computeSpringConstant(const Framework::Node* const firstNode, 
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


//////////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithmFRQ2::computeSpringConstantCenter(const Framework::Node* const firstNode, 
						     const Framework::Node* const secondNode) 
  {
  CFAUTOTRACE;
  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  typedef CFMultiMap<CFuint, CFuint> MapNodeCell;
   
  typedef MapNodeCell::MapIterator mapNodeCellIt;
  typedef MapNodeCell::MapIterator mapNodeCellItN;

  bool found = false;
  std::pair<mapNodeCellIt,mapNodeCellIt > ite=m_mapNodeCell1.find(firstNode->getLocalID(), found);
  cf_assert(found);

   
  bool foundN = false;
  std::pair<mapNodeCellItN,mapNodeCellItN > iteN=m_mapNodeCell1.find(secondNode->getLocalID(), foundN);
  cf_assert(foundN);

  CFuint cellInCommonID;

  //////////cout << "NEW COMPUTATION" << endl;
 
  /*////////cout << "FIRST NODE " << (*firstNode) << endl;
  ////////cout << "SECOND NODE " << (*secondNode) << endl;*/
  

  for (mapNodeCellIt itE = ite.first; itE != ite.second; ++itE) {
	      	   for (mapNodeCellItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
			if (itE->second == itNe->second){
				cellInCommonID = itNe->second;
				//////////cout << "CELL ID " << cellInCommonID << endl;
			}
		    }
  }
  
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
  Framework::MeshDataStack::getActive()->getTrs("InnerCells");

  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);

  geoData.idx = cellInCommonID;
  m_cell = m_cellBuilder->buildGE();

  m_cellStates = m_cell->getStates();

  for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	m_firstState[iVar] = 0;
 	m_secondState[iVar] = 0;
 }

  //s////////cout << "NEW NODE " << endl;

  for (CFuint iNode = 0; iNode <m_vecNodeCoords.size(); iNode++){
	//////////cout << "COORD NODES " << m_intCell->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) << endl;
	if (m_cell->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*firstNode)){
		//////////cout << "GG 1 " << m_vecNodeCoords[iNode] << endl;
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
			for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	      			m_firstState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
			}      
		 }

	}
	if (m_cell->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*secondNode)){
		//////////cout << "GG 2 " << m_vecNodeCoords[iNode] << endl;
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
			for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      			m_secondState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
			}      
		 }
	}
	else{
		//////////cout << "NUL" << endl;
	}

  }
  m_cellBuilder->releaseGE();
  
  CFreal springConstant = std::abs(m_secondState[m_monitorVarID] - m_firstState[m_monitorVarID]); 

  /*if ((m_order%2) == 0){
	springConstant = ((m_order+2)/2)*std::abs(m_secondState[m_monitorVarID] - m_firstState[m_monitorVarID]); 
  }
  else if ((m_order%2) != 0){
	springConstant = ((m_order+3)/2)*std::abs(m_secondState[m_monitorVarID] - m_firstState[m_monitorVarID]); 
  }*/
  
  //////////cout << "END CENTER" << endl;
  return springConstant;
}
///////////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithmFRQ2::computeSpringConstantInnerFace(const Framework::Node* const firstNode, 
						     const Framework::Node* const secondNode) 
  {
  CFAUTOTRACE;
  //////////cout << "OK" << endl;
  //////////cout << "firstNode" << firstNode->getLocalID() << endl;
  //////////cout << "secondNode" << secondNode->getLocalID() << endl;

  CFuint faceInCommonID = getFaceInCommonID(firstNode, secondNode); //get the face linking the two nodes together

  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();

  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");


  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder1->getDataGE();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
for (CFuint iNode= 0; iNode<nodes.size();iNode++){
  ////cout<<"nodes"<<nodes.size()<<endl;

}
  geoDataFace.cellsTRS = cells;
    ////cout<<"here 1"<<endl;

  geoDataFace.facesTRS = faces;
    ////cout<<"here 11"<<endl;

  geoDataFace.isBoundary = false;
    ////cout<<"here 111"<<endl;

  geoDataFace.idx = faceInCommonID;
    ////cout<<"here 1111"<<endl;

  m_face = m_faceBuilder1->buildGE();
  ////cout<<"here 1"<<endl;

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
 // reset states in nodal points
 for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	m_firstLeftState[iVar] = 0;
 	m_secondLeftState[iVar] = 0;
 	m_firstRightState[iVar] = 0;
 	m_secondRightState[iVar] = 0;
 }

    //////////cout << "NEW NODE" << endl;
    for (CFuint iNode = 0; iNode <m_vecNodeCoords.size(); iNode++){
	if (m_cells[LEFT ]->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*firstNode)){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
			for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	      			m_firstLeftState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[LEFT])[iSol])[iVar];
	                  }      
	          }

	}
	if (m_cells[LEFT]->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*secondNode)){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      			m_secondLeftState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[LEFT])[iSol])[iVar];
	                  }    
  
	          }
	}
	if (m_cells[RIGHT]->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*firstNode)){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      			m_firstRightState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[RIGHT])[iSol])[iVar];
	                  }      
	          }
	}
	if (m_cells[RIGHT]->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*secondNode)){
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
	         	for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      			m_secondRightState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_states[RIGHT])[iSol])[iVar];
	                  }      
	          }
	}
    }

  RealVector FirstNodeFluxVector(m_dim);
  RealVector SecondNodeFluxVector(m_dim);
  
  RealVector NodalVector(m_dim);
  RealVector FluxVector(m_dim);

  CFuint nbrInsideFlxPnts = 0;

  for (CFuint iDim=0 ; iDim<m_dim; iDim++){
	NodalVector[iDim] = (*secondNode)[iDim] - (*firstNode)[iDim];
  }

  //////////cout << "xfirst " << (*firstNode)[0] << "yfirst " << (*firstNode)[1] << endl;
 //////////cout << "xsecond " << (*secondNode)[0] << "ysecond " << (*secondNode)[1] << endl;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
      m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
      DistanceFirstNodeFlux[iFlx] = 1000000; 
	//////////cout << "state FLX PNT " << *(m_cellStatesFlxPnt[LEFT][iFlx]) << endl;
  }


  //////////cout << "Nodal distance " << NodalVector.norm2() << endl;

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
  	for (CFuint iDim=0 ; iDim<m_dim; iDim++){
		FirstNodeFluxVector[iDim] = m_flxPntCoords[iFlx][iDim] - (*firstNode)[iDim];
	 	SecondNodeFluxVector[iDim] = m_flxPntCoords[iFlx][iDim] - (*secondNode)[iDim];
	}
	//////////cout << "Second Node Flux Vector distance " << SecondNodeFluxVector.norm2() << endl;
	if ((FirstNodeFluxVector.norm2() < NodalVector.norm2()) && (SecondNodeFluxVector.norm2() < NodalVector.norm2())){
		DistanceFirstNodeFlux[iFlx] = FirstNodeFluxVector.norm2(); 
		nbrInsideFlxPnts +=1;
		//////////cout << "DIST" << FirstNodeFluxVector.norm2() << endl;
	}
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    *(m_actualStatesFlxPnt[LEFT][iFlx]) = 0;
    *(m_actualStatesFlxPnt[RIGHT][iFlx]) = 0;
  }


  CFreal minimumDistance;
  CFuint minimumIndex;
  std::vector<CFuint > usedIndex;
  bool alreadyUsed;
  //////////cout << "nbgInside " << nbrInsideFlxPnts << endl;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
      //////////cout << "dist first node flux " << DistanceFirstNodeFlux[iFlx] << endl;
  }

  for (CFuint insideFlx = 0; insideFlx < nbrInsideFlxPnts; ++insideFlx){
  	minimumDistance = 10000;
	for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
		alreadyUsed = false;
		for (CFuint i = 0; i<usedIndex.size(); ++i){
			if (usedIndex[i] == iFlx){
				alreadyUsed = true;
			}
		}
		if ((DistanceFirstNodeFlux[iFlx] < minimumDistance) && (!alreadyUsed)){
			minimumDistance = DistanceFirstNodeFlux[iFlx];
			minimumIndex = iFlx;
		}
	}
	usedIndex.push_back(minimumIndex);
	*(m_actualStatesFlxPnt[LEFT][insideFlx]) = *(m_cellStatesFlxPnt[LEFT][minimumIndex]);
	*(m_actualStatesFlxPnt[RIGHT][insideFlx]) = *(m_cellStatesFlxPnt[RIGHT][minimumIndex]);
  }

  CFreal stiffFirstFlux = std::abs((*(m_actualStatesFlxPnt[LEFT][0]))[m_monitorVarID] - m_firstLeftState[m_monitorVarID]) + 
                           std::abs((*(m_actualStatesFlxPnt[RIGHT][0]))[m_monitorVarID] - m_firstRightState[m_monitorVarID]);


  // compute the stiffnesses between all the flux pointss
  CFreal stiffInside = 0.;
  for (CFuint iFlxPnt = 0; iFlxPnt < nbrInsideFlxPnts-1; ++iFlxPnt){ 
	stiffInside += std::abs((*(m_actualStatesFlxPnt[LEFT][iFlxPnt]))[m_monitorVarID] - (*(m_actualStatesFlxPnt[LEFT][iFlxPnt+1]))[m_monitorVarID]) + 
                        std::abs((*(m_actualStatesFlxPnt[RIGHT][iFlxPnt]))[m_monitorVarID] - (*(m_actualStatesFlxPnt[RIGHT][iFlxPnt+1]))[m_monitorVarID]);
  }
  // compute the stiffness between the second node and the last flux point
  CFreal stiffFluxSecond = std::abs((*(m_actualStatesFlxPnt[LEFT][nbrInsideFlxPnts-1]))[m_monitorVarID] - m_secondLeftState[m_monitorVarID]) + 
                           std::abs((*(m_actualStatesFlxPnt[RIGHT][nbrInsideFlxPnts-1]))[m_monitorVarID] - m_secondRightState[m_monitorVarID]);


  CFreal stiffness = pow(stiffFirstFlux + stiffInside + stiffFluxSecond, 1.);
  m_faceBuilder1->releaseGE();


////cout<<"First node " << firstNode->getLocalID()<<endl;
////cout<<"Second node " << secondNode->getLocalID()<<endl;
////cout<<"Stiffness "<< stiffness << endl;
////cout<<" xF "<< (*firstNode)[XX]<<endl;
////cout<<" yF "<< (*firstNode)[YY]<<endl;
////cout<<" xS "<< (*secondNode)[XX]<<endl;
////cout<<" yS "<< (*secondNode)[YY]<<endl;
////cout<<"faceInCommonID "<<faceInCommonID<<endl;

////cout<<"-------------------"<<endl;
  return stiffness;
}


//////////////////////////////////////////////////////////////////////////////////
CFreal MeshFittingAlgorithmFRQ2::computeSpringConstantBoundaryFace(const Framework::Node* const firstNode, 
						     const Framework::Node* const secondNode) 
  {

  CFAUTOTRACE;

  CFuint faceInCommonID;
  CFuint commonTRS;
  
  typedef CFMultiMap<CFuint, CFuint> MapNodeFaceTRS;
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  typedef CFMultiMap<CFuint, CFuint> MapNodeTRS;

  typedef MapNodeNode::MapIterator mapItNode;
   
  typedef MapNodeFaceTRS::MapIterator mapIt;
  typedef MapNodeFaceTRS::MapIterator mapItN;

  typedef MapNodeTRS::MapIterator mapItTRS;
  typedef MapNodeTRS::MapIterator mapItNTRS;

  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  bool found = false;
  std::pair<mapIt,mapIt > ite=m_mapNodeFaceTRS1.find(firstNode->getLocalID(), found);
  cf_assert(found);

  CFreal firstNodeValue = 0;
  CFreal secondNodeValue = 0;


  bool foundN = false;
  std::pair<mapItN,mapItN > iteN=m_mapNodeFaceTRS1.find(secondNode->getLocalID(), foundN);
  cf_assert(foundN);

  for (mapIt itE = ite.first; itE != ite.second; ++itE) {
	      	   for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
			if (itE->second == itNe->second){
				faceInCommonID = itNe->second;
				//////////cout << faceInCommonID << endl;
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
  //////////cout << "NBTRS" << nbTRs << endl;
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
      //////////cout << *(m_cellStatesFlxPntBnd[iFlxPnt]) << endl;
    }
  }
 }

 //////////cout << (*m_cellStates)[1]->getLocalID()  << endl;

 std::vector< Framework::Node*  >* m_Nodes = m_intCell->getNodes();
 //////////cout << "State" << ((*(m_cellStatesFlxPntBnd[0]))[0])<< endl;
 ///////////cout << "NODES size" << m_Nodes->size() << endl;

 for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	m_firstState[iVar] = 0;
 	m_secondState[iVar] = 0;
 }
  
 //////////cout << "NEW NODE" << endl;

 for (CFuint iNode = 0; iNode <m_vecNodeCoords.size(); iNode++){
	if (m_intCell->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*firstNode)){
		//////////cout << "GG 1" << endl;
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
			for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar){
	      			m_firstState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
			}      
		 }

	}
	if (m_intCell->computeCoordFromMappedCoord(m_vecNodeCoords[iNode]) == (*secondNode)){
		//////////cout << "GG 2" << endl;
		for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol){
			for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	      			m_secondState[iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
			}      
		 }
	}
	else{
		//////////cout << "NUL" << endl;
	}

 }

  RealVector FirstNodeFluxVector(m_dim);
  RealVector SecondNodeFluxVector(m_dim);
  
  RealVector NodalVector(m_dim);
  RealVector FluxVector(m_dim);

  CFuint nbrInsideFlxPnts = 0;

  for (CFuint iDim=0 ; iDim<m_dim; iDim++){
	NodalVector[iDim] = (*secondNode)[iDim] - (*firstNode)[iDim];
  }

  //////////cout << "xfirst " << (*firstNode)[0] << "yfirst " << (*firstNode)[1] << endl;
 //////////cout << "xsecond " << (*secondNode)[0] << "ysecond " << (*secondNode)[1] << endl;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
      m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
      DistanceFirstNodeFlux[iFlx] = 1000000; 
 	//////////cout << "COORDS FLX PNT " << m_flxPntCoords[iFlx] << endl;
	//////////cout << "state FLX PNT " << *(m_cellStatesFlxPnt[LEFT][iFlx]) << endl;
  }



  //////////cout << "Nodal distance " << NodalVector.norm2() << endl;

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
  	for (CFuint iDim=0 ; iDim<m_dim; iDim++){
		FirstNodeFluxVector[iDim] = m_flxPntCoords[iFlx][iDim] - (*firstNode)[iDim];
	 	SecondNodeFluxVector[iDim] = m_flxPntCoords[iFlx][iDim] - (*secondNode)[iDim];
	}
	//////////cout << "First Node Flux Vector distance " << FirstNodeFluxVector.norm2() << endl;
	//////////cout << "Second Node Flux Vector distance " << SecondNodeFluxVector.norm2() << endl;
	if ((FirstNodeFluxVector.norm2() < NodalVector.norm2()) && (SecondNodeFluxVector.norm2() < NodalVector.norm2())){
		DistanceFirstNodeFlux[iFlx] = FirstNodeFluxVector.norm2(); 
		nbrInsideFlxPnts +=1;
		//////////cout << "DIST" << FirstNodeFluxVector.norm2() << endl;
	}
	
  }

  //////////cout << "NBR INSIDE" << nbrInsideFlxPnts << endl;

  CFreal minimumDistance;
  CFuint minimumIndex;
  std::vector<CFuint > usedIndex;
  bool alreadyUsed;

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    *(m_actualStatesFlxPntBnd[iFlx]) = 0;
  }


  for (CFuint insideFlx = 0; insideFlx < nbrInsideFlxPnts; ++insideFlx){
  	minimumDistance = 10000;
	for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx){
		alreadyUsed = false;
		for (CFuint i = 0; i<usedIndex.size(); ++i){
			if (usedIndex[i] == iFlx){
				alreadyUsed = true;
			}
		}
		if ((DistanceFirstNodeFlux[iFlx] < minimumDistance) && (!alreadyUsed)){
			minimumDistance = DistanceFirstNodeFlux[iFlx];
			minimumIndex = iFlx;
		}
	}
	usedIndex.push_back(minimumIndex);
	*(m_actualStatesFlxPntBnd[insideFlx]) = *(m_cellStatesFlxPntBnd[minimumIndex]);
  }
  
  CFreal stiffFirstFlux = std::abs((*(m_actualStatesFlxPntBnd[0]))[m_monitorVarID] - m_firstState[m_monitorVarID]);

  CFreal stiffInside = 0;
  for (CFuint iFlxPnt = 0; iFlxPnt < nbrInsideFlxPnts-1; ++iFlxPnt){ 
	stiffInside += std::abs((*(m_actualStatesFlxPntBnd[iFlxPnt]))[m_monitorVarID] - (*(m_actualStatesFlxPntBnd[iFlxPnt+1]))[m_monitorVarID]);
  }
  
  CFreal stiffFluxSecond = std::abs((*(m_actualStatesFlxPntBnd[nbrInsideFlxPnts-1]))[m_monitorVarID] - m_secondState[m_monitorVarID]);

  CFreal stiffness = 2.*pow((stiffFirstFlux + stiffInside + stiffFluxSecond), 1.);



  m_faceBuilder->releaseGE();

  return stiffness;
}




//////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithmFRQ2::truncateSpringConstant(const CFreal springConstant){
  const CFreal maxLimit = m_springTruncationData.maxLimit;
  const CFreal minLimit = m_springTruncationData.minLimit;
  const CFreal mean     = m_springTruncationData.mean;

  const CFreal truncatedSpringConstant = std::max(std::min(maxLimit/mean, springConstant/mean), minLimit/mean );
  return truncatedSpringConstant;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::solveLinearSystem(){
  assembleLinearSystem();
  m_lss->solveSys();
	
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::resizeSystemSolverToNodalData(){
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint rhsSize = std::max(nodes.size()*totalNbEqs,states.size()*totalNbEqs);
  if (rhs.size()!=rhsSize) rhs.resize(rhsSize);
  rhs = 0.;
  jacobMatrix->resetToZeroEntries();
  //////////cout<<" -----------------state size " << states.size()<< endl;
  //////////cout<<" -----------------RHS size " << rhsSize<< endl;
  //////////cout<<" -----------------totalNbEqs size " << totalNbEqs<< endl;
  //////////cout<<" -----------------nodes size " << nodes.size()<< endl;

} 
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::resizeSystemSolverToStateData(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  if (rhs.size()/totalNbEqs != states.size()) rhs.resize(states.size()*totalNbEqs); 
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::computeNodeStates()
{

  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector counter(nodalStates.size()); counter =0.; // Need the number of attached cells to a node !!
  for (CFuint i=0 ; i<nodalStates.size() ; ++i){
	  nodalStates[i] = 0.;
  }

////////////cout<<" nodal states size()  "<< nodalStates[10].size() << endl;
////////////cout<<"  m_solPolyValsAtNodes "<<m_solPolyValsAtNodes[10].size()<<endl;
  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    
    std::vector< Framework::Node*  >* m_cellNodes = currCell->getNodes();
    const CFuint nbNodes = m_cellNodes->size(); 
    ////////////cout<<"   nbNodes   "<< nbNodes << endl;
    //const CFuint nbNodes = m_cellNodes->size(); 

    //const CFuint nbNodes =  nodes.size();  //m_cellNodes->size(); 
    // get the states in this cell
    std::vector< Framework::State* >* m_cellStates = currCell->getStates();

    for (CFuint iNode = 0; iNode < nbNodes; ++iNode)  // m_nbrNodesElem = 4 for P1 BUT I want the solution evaluated at the geometric nodes
      {        
	// extrapolate the left and right states to the flx pnts
	for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	  {
      //counter[(nodes)[iNode]->getLocalID()]+=1.;
	    for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
       
	      nodalStates[(*m_cellNodes)[iNode]->getLocalID()][iVar] += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];

	      (*(m_nodalCellStates[iCell][(*m_cellNodes)[iNode]->getLocalID()]))[iVar]  += (m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
	    }      
	  }
      
      }
    m_cellBuilder->releaseGE();
    
  }

/*for(CFuint iNode=0; iNode<nodes.size(); ++iNode){
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  typedef MapNodeNode::MapIterator mapIt;
  bool found = false;
  std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find((nodes)[iNode]->getLocalID(), found);
  cf_assert(found);
  for (mapIt it = ite.first; it != ite.second; ++it) {
    counter[(nodes)[iNode]->getLocalID()]+=1.;
    //////////cout<<" counterICI " << counter[(nodes)[iNode]->getLocalID()]<< endl;
  }
}*/
  for (CFuint iNode=0 ; iNode<nodalStates.size() ; ++iNode){

    nodalStates[(nodes)[iNode]->getLocalID()] = nodalStates[(nodes)[iNode]->getLocalID()]/(m_nbOfNeighborCellsToaNode[(nodes)[iNode]->getLocalID()]);//counter[(nodes)[iNode]->getLocalID()];
  }
}
///////////////////////////////////////////////////////////////////

RealVector  MeshFittingAlgorithmFRQ2::computeIntersection(const  Framework::Node* const  a,const  Framework::Node* const  c,
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
  
CFreal MeshFittingAlgorithmFRQ2::computeElementArea2dQuads(const  Framework::Node* const  a,const  Framework::Node* const  c,
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

CFreal MeshFittingAlgorithmFRQ2::computeConstantquads(const  RealVector xyi , 
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
  //////////cout << "nbDims" << nbDims << endl;
  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    vector1[iDim]= (xyi)[iDim]-(*firstNode)[XX+iDim];
    vector2[iDim]= (xyi)[iDim]-(*secondNode)[XX+iDim];
    
    vector3[iDim]= (*thirdNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector4[iDim]= (*thirdNode)[XX+iDim]-(*secondNode)[XX+iDim];

    vector5[iDim]= (*fourthNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector6[iDim]= (*fourthNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }
  //////////cout << "oui" << endl;
  MathTools::MathFunctions::crossProd(vector1, vector2,res);
  const CFreal  area = (res.norm2()/2);
  //////////cout << "non" << endl;
  MathTools::MathFunctions::crossProd(vector3, vector4,res);
  const CFreal  area1 = (res.norm2()/2);

  MathTools::MathFunctions::crossProd(vector5, vector6,res);
  const CFreal  area2 = (res.norm2()/2);


  CFreal k = (area/elementArea)*(vector1.norm2()*vector2.norm2()*vector1.norm2()*vector2.norm2())/(4*area*area);

  CFreal k1 = (area1/elementArea)*(vector3.norm2()*vector3.norm2()*vector4.norm2()*vector4.norm2())/(4*area1*area1);

  CFreal k2 =(area2/elementArea)*(vector5.norm2()*vector5.norm2()*vector6.norm2()*vector6.norm2())/(4*area2*area2);
  // if( m_thetaMid ){
    // choose this option for a semi torsional spring analogy based on the middle angle only
  torsionConstant = k;
  // }
  //else{
    // choose this option to stiffer the mesh or on a pave mesh
    //torsionConstant =((std::max(std::max(k1,k2),k)));
    // }
  return torsionConstant;
}

//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::assembleLockedNode(const Framework::Node* node){
     Framework::DataHandle< CFreal > stiffness = socket_stiffness.getDataHandle();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;
   stiffness[node->getLocalID()]=0.;
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, 1.);
  }
  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = (*node)[XX+iDim];
  }
}

//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::assembleinRegionNode2DQuads(const Framework::Node* node){
   Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  //////////cout << "NodeDist" << nodeDistance[node->getLocalID()] << endl;
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
  //////////cout << "NodeISAD" << nodeisAD[node->getLocalID()] << endl;
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
  CFuint sumN = 0;
  for (mapIt it = ite.first; it != ite.second; ++it) {
	sumN +=1;
  }
  //////////cout << "sumN" << sumN << endl;
  cf_assert(found);
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
  //////////cout << "neigboring size" << neighboringNodes.size() << endl;
  
  std::vector<Framework::Node*>::const_iterator it;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
   	const  Framework::Node* neighborNode = *it;
	//////////cout << "xSEC" << (*neighborNode)[0] << "ySEC " << (*neighborNode)[1] << endl;
  }
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator itN1;
  std::vector<Framework::Node*>::const_iterator itSa;
  //////////cout << "xfirst " << (*node)[0] << "yfirst " << (*node)[1] << endl;
  //////////cout << "xsec " << (*secondNode)[0] << "ysec" << (*secondNode)[1] << endl;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){ // loop on the nodes 
   for (mapIt itE = ite.first; itE != ite.second; ++itE) {
    if((*it)->getLocalID() == itE->second){
    CFreal torsionConstant=0;
    std::vector<const Framework::Node*> sharedNodes;
    sharedNodes.clear();
    const  Framework::Node* neighborNode = *it;
    //////////cout << "xsecond " << (*neighborNode)[0] << "y second " << (*neighborNode)[1] << endl;
    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
	const Framework::Node* neighborNodeOfN = *itN;
         //////////cout << "xSECNNN" << (*neighborNodeOfN )[0] << "ySECNNN " << (*neighborNodeOfN )[1] << endl;
    }
    //////////cout << "neigboringN size" << neighboringNodesOfN.size() << endl;
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){ // loop on the neighbours of the neighbour node
      const Framework::Node* neighborNodeOfN = *itN;
      
      const std::vector<Framework::Node*>& neighboringNodes2 = m_edgeGraph.getNeighborNodesOfNode(node);
       for(itN1=neighboringNodes2.begin(); itN1 != neighboringNodes2.end(); ++itN1){
	const Framework::Node* neighborNodeNewloop = *itN1;
	if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
	  sharedNodes.push_back(neighborNodeNewloop);
           //////////cout << "x common " << (*neighborNodeNewloop)[0] << "y common " << (*neighborNodeNewloop)[1] << endl;
	  //////////cout << "sharedNodes size" << sharedNodes.size() << endl;
	}
      }
    }
    //////////cout << "SHARED NODES SIZE " << sharedNodes.size() << endl;
    if (sharedNodes.size()== 7 || sharedNodes.size() == 13){ 
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
		      if( d_node_neighbor.norm2() < d_neighbor_itSa.norm2() ){
		      	RealVector xyi=computeIntersection(neighborNode, sharedNodes[i],node , *itSa);
		      	CFreal elementArea=computeElementArea2dQuads(neighborNode, sharedNodes[i],node , *itSa);
		      	torsionConstant+=computeConstantquads(xyi, node , neighborNode,  sharedNodes[i]   ,   *itSa  , elementArea);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      //////////cout << "torsionConstant" << torsionConstant << endl;
      CFreal springConstant = 0;
      if (!isInSideCell(node) && !isInSideCell(neighborNode)){
	  //////////cout << "GG1" << endl;
	  springConstant =computeSpringConstantInnerFace(node,neighborNode);
	  //////////cout << "GGAFTER" << endl;
      }
      else if(isInSideCell(node) || isInSideCell(neighborNode)){
	  //////////cout << "GG2 "<< endl;
	  springConstant =computeSpringConstantCenter(node,neighborNode);
	  //////////cout << "GG2" << endl;
      }
      //////////cout << "MODIF" << endl;
      const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant); 
      //////////cout << "normalizedSpringConstant" << normalizedSpringConstant << endl;
      CFreal f = 1.;
      if (m_smoothSpringNetwork) {
      	//CFreal f = (7.-2.)/(0.02-0.01)*nodeDistance[node->getLocalID()]+ 1.2 - (7.-2.)/(0.02-0.01)*0.01;  //0.09
	// FB :test case dependent
      	 CFreal f = (5.-1.)/(0.0001-m_acceptableDistanceQ2)*nodeDistance[node->getLocalID()]+ 1. - (5.-1.)/(0.0001-m_acceptableDistanceQ2)*m_acceptableDistanceQ2;  // 5
      }
      const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
      const CFreal stiff=normalizedSpringConstant+(pow(torsionConstant,1.)*pow(normalizedSpringConstant,f));  
      //////////cout << "kstiff" << stiff << endl;
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

//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::assembleInnerNode(const Framework::Node* node){
   SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
     getTrs("InnerCells");
   const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
   // Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
   //Framework::DataHandle<Framework::Node*, Framework::GLOBAL> neighborNodes = socket_nodes.getDataHandle();
     const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
     std::vector<Framework::Node*>::const_iterator itN;

  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  typedef MapNodeNode::MapIterator mapIt;
  bool found = false;
  std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
  cf_assert(found);

   Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
   Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
   Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
   const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
   const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
   const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
   //2D quads
     CFreal sumOffDiagonalValues = 0.;
      for(itN=neighboringNodes.begin(); itN != neighboringNodes.end(); ++itN){ 
        const Framework::Node* neighborNode = *itN;
        for (mapIt it = ite.first; it != ite.second; ++it) {
          if(neighborNode->getLocalID() == it->second){
	    CFreal springConstant = 0;
	    //if (!isInSideCell(node) && !isInSideCell(neighborNode)){
		//////////cout << "GG1" << endl;
		springConstant =computeSpringConstantInnerFace(node,neighborNode);
		//////////cout << "GGAFTER" << endl;
	    //}
	    //else if(isInSideCell(node) || isInSideCell(neighborNode)){
		//////////cout << "GG2 "<< endl;
		//springConstant =computeSpringConstantCenter(node,neighborNode);
		//////////cout << "GG2" << endl;
	   // }
		 //////////cout << "spring constant" << springConstant << endl;
            const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
	           //////////cout<<" normalized spring constant  "<< normalizedSpringConstant << endl;
            sumOffDiagonalValues +=normalizedSpringConstant ;
            const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
	          const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
            //////////cout<< "IN: Global node ID RG  " << rowGlobalID << "   Local ID    " << node->getLocalID()<<endl;
            //////////cout<< "IN: Global node ID CG  " << colGlobalID << "   Local ID    " << neighborNode->getLocalID()<<endl;
            //////////cout<< "idxMapping size        "<< idxMapping.getRowID(node->getLocalID()) << endl;
            //////////cout<< "IN: Global node ID CG  " << colGlobalID << "   Local ID    " << neighborNode->getLocalID()<<endl;
	          for(CFuint iDim=0; iDim<nbDims; ++iDim){
	            jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,normalizedSpringConstant);
            }
	          stiffness[node->getLocalID()]=normalizedSpringConstant;
            CFLog(INFO,"springConstant: " << springConstant << "\n");
            CFLog(INFO,"normalizedSpringConstant: " << normalizedSpringConstant << "\n");
            cout<<"   "<<endl;

          }
        }
      }
     for(CFuint iDim=0; iDim<nbDims; ++iDim){     
       const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;
       //////////cout<< "diagvalues  "<<globalID << endl;
       jacobMatrix->addValue(globalID+iDim, globalID+iDim, -sumOffDiagonalValues);
     }

     //Right hand side
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
       //////////cout<<" RHS Filling " << endl;
       rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues);
     }
     //jacobMatrix->printToScreen();

 }
////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::assembleMovingInBoundaryNode(const Framework::Node* node){
   
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
	if (dist.norm2() < m_equilibriumSpringLength*1. && blocked == false){
		////////cout << "blocked" << endl;
		blocked = false;
	}
    }
  } 
  
  if (!blocked){
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
	//////////cout << "Moving" << endl;
	CFreal sumOffDiagonalValues = 0.;
  	typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
     	typedef MapNodeNode::MapIterator mapIt;
     	bool found = false;
     	std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
     	cf_assert(found);  
	const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
         std::vector<Framework::Node*>::const_iterator itN;
  	for(itN=neighboringNodes.begin(); itN != neighboringNodes.end(); ++itN){
	   //////////cout << "neighbor" << endl;
	   Framework::Node* neighborNode = *itN;
	   for (mapIt it = ite.first; it != ite.second; ++it) {
               if (neighborNode->getLocalID() == it->second){
		const bool neighborIsBoundary = m_boundaryNodes.find(neighborNode) != m_boundaryNodes.end();
		CFreal springConstant = 0;
		if (neighborIsBoundary){
			      CFuint commonTRS;
            bool sameTRS = false;

            typedef CFMultiMap<CFuint, CFuint> MapNodeTRS1;
            typedef MapNodeTRS1::MapIterator mapItTRS;
            typedef MapNodeTRS1::MapIterator mapItNTRS;
            bool foundTRS = false;
            std::pair<mapItTRS,mapItTRS > itTRS =m_mapNodeTRS1.find(node->getLocalID(), foundTRS);
            cf_assert(foundTRS);
            bool foundNTRS = false;
            std::pair<mapItNTRS,mapItNTRS > itNTRS =m_mapNodeTRS1.find(neighborNode->getLocalID(), foundNTRS);
            cf_assert(foundNTRS);
            for (mapItTRS ittrs = itTRS.first; ittrs != itTRS.second; ++ittrs) {
	            for (mapItNTRS itntrs = itNTRS.first; itntrs != itNTRS.second; ++itntrs) {
		            if (ittrs->second == itntrs->second){
		              commonTRS = itntrs->second;
                  sameTRS = true;
			          }
		          }
            }
            if (sameTRS){
				      springConstant =computeSpringConstantBoundaryFace(node, neighborNode);
            }
            else{
              springConstant =computeSpringConstantInnerFace(node, neighborNode);
            }
			//////////cout << "spring constant" << springConstant << endl;
		}
		//else if (isInSideCell(neighborNode)){
			//////////cout << "NO" << endl;
			//////////cout << "xfirst " << (*node)[0] << "yfirst " << (*node)[1] << endl;
			//////////cout << "xsec " << (*neighborNode)[0] << "ysec" << (*neighborNode)[1] << endl;
    	//		springConstant = computeSpringConstantCenter(node, neighborNode); 
		//}
		else {
			springConstant = computeSpringConstantInnerFace(node, neighborNode);
		}
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
    		//CFreal normalizedSpringConstant =  pow(physicalSpringConstant,m_order);
		CFreal normalizedSpringConstant =  physicalSpringConstant;
		//////////cout << "it is normalized" << endl;
		RealVector dist(nbDims);
                  for (CFuint iDim=0; iDim<nbDims; ++iDim){
		dist[iDim] = (*neighborNode)[XX+iDim]-(*node)[XX+iDim];
		}
		const CFreal stiff= normalizedSpringConstant;
          	stiffness[node->getLocalID()]=stiff;
		sumOffDiagonalValues += stiff;
		//////////cout << "it is stiff" << endl;
    		const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
    		const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
    		for(CFuint iFreeDim=0; iFreeDim<freeDims.size(); ++iFreeDim){
      			jacobMatrix->addValue(rowGlobalID+freeDims[iFreeDim], colGlobalID+freeDims[iFreeDim], stiff);
    		}
		//////////cout << "it is jacob" << endl;
	   }
          }
  	}
	//////////cout << "end of this node" << endl;
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
void MeshFittingAlgorithmFRQ2::assembleLinearSystem(){
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
 //////////cout << "NB nodes " << nbNodes << endl;
 for (CFuint iNode = 0; iNode < nbNodes; ++iNode){
         //////////cout << "iNode" << iNode << endl;
	if (!nodes[iNode]->isParUpdatable()){ 
	  //do nothing
	}
  else{

	if(insideRegion(nodes[iNode]) ){
              // assembleInnerNode(nodes[iNode]);
	     //assembleinRegionNode2DQuads(nodes[iNode]);
		assembleLockedNode(nodes[iNode]);

	////////cout<<" m_acceptableDistanceQ2   "<< m_acceptableDistanceQ2 << endl;
	}
	else{ 
	if((isNodeMovingInBoundary(nodes[iNode])) ){
	  //////////cout << "moving" << endl;	
	  //////////cout << "xfirst " << (*nodes[iNode])[0] << "yfirst " << (*nodes[iNode])[1] << endl;
	  assembleMovingInBoundaryNode(nodes[iNode]);
	  //assembleLockedNode(nodes[iNode]);
	}

	else{
	    if(isNodeLocked(nodes[iNode])){
	      assembleLockedNode(nodes[iNode]);
	    }
	    else{
	      //////////cout << "inner" << endl;	
	      //////////cout << "xfirst " << (*nodes[iNode])[0] << "yfirst " << (*nodes[iNode])[1] << endl;
	      //////////cout << "before inner" << endl;
	      assembleInnerNode(nodes[iNode]);
	      //////////cout << "after inner" << endl;
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
CFreal MeshFittingAlgorithmFRQ2::determinantOfMatrix(std::vector< RealVector >  mat) 
{ 
    CFreal ans; 
    ans = mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1] * mat[1][2]) - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]); 
	 CFLog(NOTICE," ans  " << ans<< "\n");
    return ans; 
} 
/////////////////////////////////////////////////////////////////////////////
RealVector MeshFittingAlgorithmFRQ2::findSolution(std::vector< RealVector > coeff) 
{ 

	RealVector a_b_c (3); a_b_c=0.;
    // Matrix d using coeff as given in cramer's rule 
	vector< RealVector > d;  d.resize(3);
	vector< RealVector > d1; d1.resize(3);
	vector< RealVector > d2; d2.resize(3);
	vector< RealVector > d3; d3.resize(3);
	for (CFuint i=0; i<3; ++i){
		d[i].resize(3);
		d1[i].resize(3);
		d2[i].resize(3);
		d3[i].resize(3);
	}
/*    CFreal d[3][3] = { 
        { coeff[0][0], coeff[0][1], coeff[0][2] }, 
        { coeff[1][0], coeff[1][1], coeff[1][2] }, 
        { coeff[2][0], coeff[2][1], coeff[2][2] }, 
    }; */

	for(CFuint i=0; i<3 ; ++i){
		for(CFuint j=0; j<3; ++j){
			d[i][j]=coeff[i][j];
		}
	}
    // Matrix d1 using coeff as given in cramer's rule 
   /* CFreal d1[3][3] = { 
        { coeff[0][3], coeff[0][1], coeff[0][2] }, 
        { coeff[1][3], coeff[1][1], coeff[1][2] }, 
        { coeff[2][3], coeff[2][1], coeff[2][2] }, 
    }; */


	for( CFuint i=0; i<3; ++i){
		d1[i][0] = coeff[i][3];
		d1[i][1] = coeff[i][1];
		d1[i][2] = coeff[i][2];
	}
	

  for( CFuint i=0; i<3; ++i){
		d2[i][0] = coeff[i][0];
		d2[i][1] = coeff[i][3];
		d2[i][2] = coeff[i][2];
	}


  for( CFuint i=0; i<3; ++i){
		d3[i][0] = coeff[i][0];
		d3[i][1] = coeff[i][1];
		d3[i][2] = coeff[i][3];
	}

/*

    // Matrix d2 using coeff as given in cramer's rule 
    CFreal d2[3][3] = { 
        { coeff[0][0], coeff[0][3], coeff[0][2] }, 
        { coeff[1][0], coeff[1][3], coeff[1][2] }, 
        { coeff[2][0], coeff[2][3], coeff[2][2] }, 
    }; 
    // Matrix d3 using coeff as given in cramer's rule 
    CFreal d3[3][3] = { 
        { coeff[0][0], coeff[0][1], coeff[0][3] }, 
        { coeff[1][0], coeff[1][1], coeff[1][3] }, 
        { coeff[2][0], coeff[2][1], coeff[2][3] }, 
    }; */
  // Calculating Determinant of Matrices d, d1, d2, d3 

    CFreal D  = determinantOfMatrix(d); 
    CFreal D1 = determinantOfMatrix(d1); 
    CFreal D2 = determinantOfMatrix(d2); 
    CFreal D3 = determinantOfMatrix(d3); 

    // Case 1 
    if (D != 0) { 
        // Coeff have a unique solution. Apply Cramer's Rule 
        CFreal a = D1 / D; 
        CFreal b = D2 / D; 
        CFreal c = D3 / D; // calculating z using cramer's rule 
        a_b_c[0]=a;
				a_b_c[1]=b;
				a_b_c[2]=c;
	 CFLog(NOTICE," yes ok   " << "\n");
    } 
    // Case 2 
    else{
	 CFLog(NOTICE," infiniteSolution  " << "\n");
		}

  return  a_b_c;
} 
/////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::updateNodePositions () {
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::updateNodePositions()\n"); 
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
  Framework::DataHandle < CFreal > rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  vector< RealVector > halfStep; halfStep.resize(nodes.size());
  vector < CFuint > boundaryNodeID; boundaryNodeID.resize(0);

 /* for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
	//if (nodes[iNode]->isParUpdatable()){ 
	  
    if (nodes[iNode]->isParUpdatable()) {
      Framework::Node& currNode = *nodes[iNode];
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
        currNode[XX+iDim] = currNode[XX+iDim]*(1.-m_meshAcceleration) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration; 
      }
}
    //}
  }*/


  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
   CFreal f=1.;
    if (nodes[iNode]->isParUpdatable()) {
      Framework::Node& currNode = *nodes[iNode];
      bool exit = false;
  /*    if (m_smoothNodalDisp){
	// FB : test case dependent

	if( nodeDistance[nodes[iNode]->getLocalID()] < 1. &&  insideRegion(nodes[iNode])==false ){ 				
	  f=(1./(1.-m_acceptableDistanceQ2)) * nodeDistance[nodes[iNode]->getLocalID()] - (m_acceptableDistanceQ2/(1.-m_acceptableDistanceQ2)); 
	}
      }*/
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
	currNode[XX+iDim] =  currNode[XX+iDim]*(1.-m_meshAcceleration*f) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration*f;
      }
    }
  }


/// Middle point treatment

/*

  for(CFuint iNode =0; iNode < nodes.size(); ++iNode){
     halfStep[nodes[iNode]->getLocalID()].resize(nbDims);
  }
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
      Framework::Node& currNode = *nodes[iNode];
	//if(nodes[iNode]->isParUpdatable()){
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
				halfStep[(nodes[iNode])->getLocalID()][XX+iDim] =  currNode[XX+iDim]*(1.-m_meshAcceleration) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration;
    //  } 
	}
    if ( isBoundaryNode(nodes[iNode])) { //nodes[iNode]->isParUpdatable() && 
		boundaryNodeID.push_back((nodes[iNode])->getLocalID());
    }
   }
  for(CFuint iNode = 0 ; iNode< nodes.size(); ++iNode){
    if (m_nbOfNeighborCellsToaNode[(nodes[iNode])->getLocalID()] == 1 && isBoundaryNode(nodes[iNode]) && nbOfConnectedFaces[(nodes[iNode])->getLocalID()] == 1){
		halfStep[(nodes[iNode])->getLocalID()]  = 0.;
		typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  		typedef MapNodeNode::MapIterator mapIt;
  		bool found = false;
  		std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find((nodes[iNode])->getLocalID(), found);
  		cf_assert(found);
        vector< RealVector > position; position.resize(3);
		for (CFuint iSize = 0; iSize< 3; ++iSize){
			position[iSize].resize(2);
		}
		CFuint i=0;
    CFreal x_centerNode;
        for (mapIt it = ite.first; it != ite.second; ++it) {
			////////cout<< "boundaryNodeID.size()   "<<ite.getSize()<<endl;
			for (CFuint iBFnode = 0; iBFnode< boundaryNodeID.size(); ++iBFnode){ //This is useful to remove from the multimap the center point
          if(boundaryNodeID[iBFnode] ==  it->second ){
					//for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
					CFreal x = halfStep[boundaryNodeID[iBFnode]][XX+0];
					CFreal y = halfStep[boundaryNodeID[iBFnode]][XX+1];
					position[i][0]=x;
					position[i][1]=y;
					i=i+1;
					//}	
				}
				else{
				x_centerNode = halfStep[it->second][XX+0];

}
		 	}
		}
		cf_assert(i==2);

		position[i][0]= halfStep[(nodes[iNode])->getLocalID()][XX+0];
		position[i][1]= halfStep[(nodes[iNode])->getLocalID()][XX+1];
	  //CFLog(NOTICE,"position[i][0]: " <<i <<"  "<< position[i][0] << "\n");
	  //CFLog(NOTICE,"position[i][1]: " <<i <<"  "<< position[i][1] << "\n");
*/
/*
System to solve to find a,b and c and get the second order equation describing the edge going through 3 boundary nodes
ax_i² + b x_i + c = y_i
where i goes from 0 to 2 
*/
/*
  // storing coefficients of linear equations in coeff matrix 
	vector< RealVector > coeff; coeff.resize(3);
	for (CFuint iSize = 0; iSize< 3; ++iSize){
		coeff[iSize].resize(4);
	}
	
   for(CFuint i =0; i< 3; ++i){
		coeff[i][0]=position[i][0]*position[i][0];
		coeff[i][1]=position[i][0];
		coeff[i][2]=1.;
		coeff[i][3]=position[i][1];
	  //CFLog(NOTICE,"coeff[i][0]  " <<i <<"  "<< coeff[i][0] << "\n");
	  //CFLog(NOTICE,"coeff[i][1]  " <<i <<"  "<< coeff[i][1] << "\n");
	  //CFLog(NOTICE,"coeff[i][2]  " <<i <<"  "<< coeff[i][2] << "\n");
	  //CFLog(NOTICE,"coeff[i][3]  " <<i <<"  "<< coeff[i][3] << "\n");
}*/

/*
	RealVector a_b_c; a_b_c.resize(3);
  a_b_c = findSolution(coeff);
	 CFLog(NOTICE," a_b_c  " << a_b_c<< "\n");
	if( a_b_c.norm2()>1e-3){
  	CFreal x_middlepoint = x_centerNode;//(position[0][0] + position[1][0])*0.5;
  	halfStep[(nodes[iNode])->getLocalID()][XX+0] = x_middlepoint ; // Take the middle of the 2 vertex points
  	CFreal y_middlepoint = a_b_c[0]* pow(x_middlepoint,2.) + a_b_c[1]* pow(x_middlepoint, 1.) +  a_b_c[2]; //evaluate the X coord within the 2nd order equation 
  	halfStep[(nodes[iNode])->getLocalID()][XX+1] = y_middlepoint;
	}
	if( a_b_c.norm2()<=1e-3){
  	CFreal x_middlepoint = (position[0][0] + position[1][0])*0.5;
  	halfStep[(nodes[iNode])->getLocalID()][XX+0] = x_middlepoint ; // Take the middle of the 2 vertex points
  	CFreal y_middlepoint = (position[0][1] + position[1][1])*0.5;
  	halfStep[(nodes[iNode])->getLocalID()][XX+1] = y_middlepoint;
	}
}
}
*/

/*
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
   CFreal f=1.;
    if (nodes[iNode]->isParUpdatable()) {
      	Framework::Node& currNode = *nodes[iNode];
      	bool exit = false;
		//if (m_smoothNodalDisp){
		if (m_nbOfNeighborCellsToaNode[(nodes[iNode])->getLocalID()] == 1 && isBoundaryNode(nodes[iNode]) && nbOfConnectedFaces[(nodes[iNode])->getLocalID()] == 1){
      		for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
				currNode[XX+iDim] = halfStep[(nodes[iNode])->getLocalID()][XX+iDim];
      		}
		}
    	else{
    //if (nodes[iNode]->isParUpdatable()) {
      		for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
				currNode[XX+iDim] =  currNode[XX+iDim]*(1.-m_meshAcceleration*f) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration*f;
			}
}
      	}
    }
   //}
*/



		
  //////////cout << "non" << endl;
  //synchronize Nodes
  nodes.beginSync();
  nodes.endSync();
}

/////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFRQ2::insideRegion(Framework::Node* node)
{
  bool inRegion=false;
  Framework::Node& currNode = *node;
  const CFuint nodeID = currNode.getLocalID();
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
 // FB : uncomment isNodeLocked(node)==false f you are not using nodal movememt interpolation
  if  ((nodeisAD[nodeID] == true)){  //&& (isNodeLocked(node)==false)){
    inRegion=true;	    
  }
  RealVector n_0 (2);
  n_0[0]= 0.00045;
  n_0[1]=0.0002;
  RealVector Node_fixedPoint (2);
  Node_fixedPoint[0]= (*node)[XX+0]-0.0001;
  Node_fixedPoint[1]= (*node)[XX+1]-0.0006;
  bool small = false;
  CFreal direction = MathTools::MathFunctions::innerProd(Node_fixedPoint, n_0);
  if (direction<0.){
	small = true;
	}
  return (inRegion);
}
//}
/////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFRQ2::isNodeLocked( Framework::Node* node)
{
  const bool isBoundary = m_boundaryNodes.find(node) != m_boundaryNodes.end();

  return (isBoundary) ;
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::findBoundaryNodes() // probably need to add already computed vector
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::createConnectivity()" << "\n");
                  typedef CFMultiMap<CFuint, CFuint> MapNodeTRS1;
                typedef MapNodeTRS1::MapIterator mapItTRS;
                typedef MapNodeTRS1::MapIterator mapItNTRS;
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
    std::vector<CFuint> alreadyComputed;
  CFuint NBNODES = 0;
  //std::vector< Framework::Node > m_boundaryNodes;
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();

  for (CFuint iCell=0; iCell<nbCells; ++iCell) {
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();

    const std::vector<Framework::GeometricEntity*  >* facesInCell = currCell->getNeighborGeos();
    Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryCell = 
      m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();

    const CFuint nbFaces = facesInCell->size();
    for (CFuint iFace=0; iFace<nbFaces; ++iFace) {
      if((*m_isFaceOnBoundaryCell)[iFace]) {
	      std::vector<Framework::Node* >* faceNodes = (*facesInCell)[iFace]->getNodes();
	      const CFuint nbFaceNodes = faceNodes->size();
	      for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {
          bool isOnBoundary = false;
  	      bool done = false;

          std::pair<mapItTRS,mapItTRS > itTRS =m_mapNodeTRS1.find((*faceNodes)[iNode]->getLocalID(), isOnBoundary);
          if (isOnBoundary){
	        for(CFuint jxx=0; jxx<alreadyComputed.size(); ++jxx){
	          if(alreadyComputed[jxx] ==  (*faceNodes)[iNode]->getLocalID()){
	            done = true;
	          }
	        }

	        if(done==false  ){

	    	    alreadyComputed.push_back((*faceNodes)[iNode]->getLocalID());
	    	    m_boundaryNodes.insert((*faceNodes)[iNode]);
		        NBNODES +=1;
	  }
	}
	      } 
	    }
    }
    m_cellBuilder->releaseGE();
  }
  //cout<<"   in find boundary nodes  "<<m_boundaryNodes.size()<<endl;
    //cout<<"  NBNODES  "<<NBNODES<<endl;

}

//////////////////////////////////////////////////////////////////////////////

bool MeshFittingAlgorithmFRQ2::isNodeMovingInBoundary(Framework::Node* node){
  
  return (m_mapNodeIDNormal.find(node->getLocalID()) != m_mapNodeIDNormal.end());
}
//////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFRQ2::isBoundaryNode(Framework::Node* node){ 
  const bool isBoundary = m_boundaryNodes.find(node) != m_boundaryNodes.end();
  return (isBoundary) ;
}


//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::findCenterNodes()
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::createConnectivity()" << "\n");
    Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();

  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  std::vector< bool > m_inFace;
  m_inFace.resize(nodes.size());
  for(CFuint iNode=0; iNode<nodes.size(); ++iNode){
    m_inFace[iNode] = false;
  }
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  //std::vector< Framework::Node > m_boundaryNodes;
  for (CFuint iCell=0; iCell<nbCells; ++iCell) {
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    const std::vector<Framework::GeometricEntity* >& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size();
    std::vector< Framework::Node*  >& m_cellNodes = *currCell->getNodes();
    const CFuint nbNodes = m_cellNodes.size(); 
    for (CFuint iFace=0; iFace<nbFaces; ++iFace) {
	    std::vector<Framework::Node* >& faceNodes = *facesInCell[iFace]->getNodes();
      for(CFuint iNodeinFace=0; iNodeinFace<faceNodes.size(); ++iNodeinFace){
        for(CFuint iNodeinCell=0; iNodeinCell<nbNodes; ++iNodeinCell){
          if(m_cellNodes[iNodeinCell] == faceNodes[iNodeinFace]){
            m_inFace[m_cellNodes[iNodeinCell]->getLocalID()]= true;
	    	    //m_CenterNodes[(*faceNodes)[iNodeinFace]->getLocalID()] = faceNodes[iNodeinFace];
          }
        }
      }
    }
  m_cellBuilder->releaseGE();
  }
  for (CFuint iNode=0; iNode<nodes.size(); ++iNode){
    if(m_inFace[nodes[iNode]->getLocalID()] == false){
        m_CenterNodes.insert(nodes[iNode]->getLocalID());
    }


  }
  ////////cout<<"   in find center nodes"<<m_CenterNodes.size()<<endl;


}


//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::nbOfConnectedFacesToaNode()
{

  //FB: Not accurate: Only for the corner node (2 Faces) and middle boundary point (1 Face), which is the goal for this function 
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::createConnectivity()" << "\n");
    Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();

  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  //std::vector< Framework::Node > m_boundaryNodes;
  for (CFuint iCell=0; iCell<nbCells; ++iCell) {
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    const std::vector<Framework::GeometricEntity* >& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size();
    std::vector< Framework::Node*  >& m_cellNodes = *currCell->getNodes();
    const CFuint nbNodes = m_cellNodes.size(); 
    for (CFuint iFace=0; iFace<nbFaces; ++iFace) {
	    std::vector<Framework::Node* >& faceNodes = *facesInCell[iFace]->getNodes();
      for(CFuint iNodeinFace=0; iNodeinFace<faceNodes.size(); ++iNodeinFace){
        nbOfConnectedFaces[faceNodes[iNodeinFace]->getLocalID()]=nbOfConnectedFaces[faceNodes[iNodeinFace]->getLocalID()]+1;
      }
    }
  m_cellBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFRQ2::isInSideCell(const Framework::Node* node){ 
  const bool isinsideCell = m_CenterNodes.find(node->getLocalID()) != m_CenterNodes.end();
  return (isinsideCell) ;
}

//////////////////////////////////////////////////////////////////////////////
   
void MeshFittingAlgorithmFRQ2::computeMovingInBoundaryNodeNormals()
{
  CFAUTOTRACE;
  const CFuint nbDim = Framework::PhysicalModelStack::getActive()->getDim();
  //Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle(); 
  std::vector<std::set<CFuint> > mapNodeID2BoundaryFaceIDs = getMapMovingInBoundaryNodeID2FaceID();
  //Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle(); 
  std::set<Framework::Node*>::iterator it;
  RealVector normal1(nbDim), normal2(nbDim);
  for(CFuint iNodeID=0; iNodeID<mapNodeID2BoundaryFaceIDs.size(); ++iNodeID){
    //Check if all normals of the connected faces are collinear up to a certain tolerance
    //that means check if the dot product of the normals are bellow that tolerance

    if(mapNodeID2BoundaryFaceIDs[iNodeID].size()>=nbDim){

      std::set<CFuint>& connectedFaceIds = mapNodeID2BoundaryFaceIDs[iNodeID];
      bool areNormalsCollinear  = true;
      std::set<CFuint>::iterator itFaceID1, itFaceID2;
      for(itFaceID1=connectedFaceIds.begin(); itFaceID1!=connectedFaceIds.end(); ++itFaceID1){
        for(itFaceID2=connectedFaceIds.begin(); itFaceID2!=connectedFaceIds.end(); ++itFaceID2){
          if(itFaceID1 != itFaceID2){
	   //////////cout << "faceDifferent" << endl;
            for(CFuint iDim=0; iDim<nbDim; ++iDim){
	      //////////cout<< "  (*itFaceID1)                "<<(*itFaceID1)<<endl;
	      //////////cout<< "  (*itFaceID1)*nbdims +idim   "<<(*itFaceID1)*nbDim+iDim<<endl;
	      //////////cout<< normals.size() << endl;
	      //////////cout << "start" << endl;
              normal1[iDim] = m_normalsAMR[(*itFaceID1)*nbDim+iDim];
              normal2[iDim] = m_normalsAMR[(*itFaceID2)*nbDim+iDim];
	     //////////cout << "normal1 " << normal1 << endl;
	     //////////cout << "normal2 " << normal2 << endl;
            }
            normal1.normalize(); normal2.normalize();
            const CFreal cosAngle = MathTools::MathFunctions::innerProd(normal1, normal2);
	    //////////cout<< " cosine angle    " << cosAngle<< endl;
            areNormalsCollinear &= (1.- std::abs(cosAngle)< 1e-2);
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
        //////////cout << "av normal" << averageNormal << endl;
        //CFLog(NOTICE,"averageNormal  "<<averageNormal<< "/n");
        m_mapNodeIDNormal.insert(addPair);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
 
std::vector<std::set<CFuint> > MeshFittingAlgorithmFRQ2::getMapMovingInBoundaryNodeID2FaceID() /// Will need to be modified to take into account curvature of elements !!
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
    //////////cout << "NBR BND Faces" << nbFaces << endl;
    for(CFuint i=0; i<nbFaces; ++i){
      facesData.idx = i;
      const Framework::GeometricEntity *const face = m_faceBuilder->buildGE();
      const std::vector<Framework::Node*> faceNodes = face->getNodes();
      const CFuint faceID = face->getID();
      const CFuint nbFaceNodes = faceNodes.size();
      //////////cout << "Nb Face Nodes" << nbFaceNodes << endl;
      for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {
	const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(faceNodes[iNode]);
	//////////cout << "neighboringSize" << neighboringNodes.size() << endl;
          const CFuint nodeID = faceNodes[iNode]->getLocalID();
          //////////cout << "xfirst " << (*faceNodes[iNode])[0] << "yfirst " << (*faceNodes[iNode])[1] << endl;
	 //nbrNodes +=1;
          (mapBoundaryNodeID2FaceIDs[nodeID]).insert(faceID);
      }
      //////////cout << "xfirst " << (*faceNodes[1])[0] << "yfirst " << (*faceNodes[1])[1] << endl;
      //////////cout << "xsec " << (*faceNodes[0])[0] << "ysec" << (*faceNodes[0])[1] << endl;
      CFreal deltaX = (*faceNodes[1])[0] - (*faceNodes[0])[0];
      CFreal deltaY = (*faceNodes[1])[1] - (*faceNodes[0])[1];
      //////////cout << "deltaX" << deltaX << endl;
      //////////cout << "deltaY" << deltaY << endl;
      m_normalsAMR[faceID*nbDim+0] = -deltaY;
      m_normalsAMR[faceID*nbDim+1] = deltaX;
      m_faceBuilder->releaseGE();
    }
  }
  //////////cout << "Nbr NODES" << nbrNodes << endl;
  return mapBoundaryNodeID2FaceIDs;

}

//////////////////////////////////////////////////////////////////////////////

 void MeshFittingAlgorithmFRQ2::createGeneralConnectivityFR()
{ 
  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR(); 
  //Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  //for(CFuint iNode=0 ; iNode<nodes.size(); ++iNode){
  //const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(nodes[iNode]);
  //std::vector<Framework::Node*>::const_iterator itN;
  //////////cout<<" neighbor node size "<< neighboringNodes.size()<< endl;
   // }
   
  m_edgeGraphN.setNodeDataSocketFR(socket_nodes);
  m_edgeGraphN.computeConnectivityFR(); 
  }


//////////////////////////////////////////////////////////////////////////////
 
void MeshFittingAlgorithmFRQ2::triggerRecomputeMeshData() {
  std::string msg;
  Common::SafePtr<Common::EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  const std::string ssname = Framework::SubSystemStatusStack::getCurrentName();   
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MESHADAPTER_AFTERMESHUPDATE"), msg );
  /*////////cout << "A0" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_UNSETUP"), msg );
  //////////cout << "A1" << endl;
 // event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_UNPLUGSOCKETS"), msg );
  //////////cout << "A2" << endl;
  // event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_DESTROYSUBSYSTEM"), msg ); 
 // ////////cout << "A3" << endl;
  // event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_BUILDSUBSYSTEM"), msg );     
 //////////cout << "A4" << endl;
  // event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_PLUGSOCKETS"), msg );
  ////////cout << "A5" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_BUILDMESHDATA"), msg );
  ////////cout << "A6" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_SETUP"), msg );
  ////////cout << "A7" << endl;*/
}
    }/// Namespace FR
  /// Namespace COOLFluiD
  
}
