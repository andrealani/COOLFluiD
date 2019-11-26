#include "Common/PE.hh"
#include "Common/EventHandler.hh"
#include "MathTools/LeastSquaresSolver.hh"
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

#include "FluxReconstructionMethod/MeshFittingAlgorithmFR.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"



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

MethodCommandProvider<MeshFittingAlgorithmFR, 
		      DataProcessingData,
		      FluxReconstructionModule>
MeshFittingAlgorithmFRFluxReconstructionProvider("MeshFittingAlgorithmFR");

///////////////////////////////////////////////////////////////////


void MeshFittingAlgorithmFR::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("minPercentile","Percentile for minimum spring value");
  options.addConfigOption< CFreal >("maxPercentile","Percentile for maximum spring value");
  options.addConfigOption< CFreal >("meshAcceleration","How fast the mesh moves in mesh steps");
  options.addConfigOption< CFuint >("monitorVarID","Monitor variable ID (from State) for mesh adaptation");
  options.addConfigOption< CFuint >("monitorPhysVarID","Monitor physical variable ID (from physical data) for mesh adaptation");
  options.addConfigOption< CFreal >("equilibriumSpringLength","Length of spring for equilibrium");
  options.addConfigOption< CFreal >("ratioBoundaryToInnerEquilibriumSpringLength","ratio between the equilibrium length of a Boundary spring to an Inner spring");
  options.addConfigOption< std::vector<std::string> >("unlockedBoundaryTRSs","TRS's to be unlocked");
}

///////////////////////////////////////////////////////////////////

MeshFittingAlgorithmFR::MeshFittingAlgorithmFR(const std::string& name) :

  Framework::DataProcessingCom(name),  
  socket_stiffness("stiffness"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_normals("normals"),
  socket_rhs("rhs"), 
  m_lss(CFNULL),
  m_geoBuilder(),
  m_cellBuilder(CFNULL),
  m_faceBuilder(),
  m_order(),
  m_frData(CFNULL),
  m_solPolyValsAtNodes(CFNULL),
  m_nbrNodesElem(),
  m_nbrSolPnts()
  
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
}

//////////////////////////////////////////////////////////////////////////////

MeshFittingAlgorithmFR::~MeshFittingAlgorithmFR()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFR::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

 Framework::DataProcessingCom::configure(args);
  
  cf_assert(m_minPercentile >= 0.);
  cf_assert(m_minPercentile < m_maxPercentile);
  cf_assert(m_maxPercentile <= 1.);
  cf_assert(m_meshAcceleration > 0. && m_meshAcceleration < 1.);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
MeshFittingAlgorithmFR::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_normals);
  result.push_back(&socket_rhs);

  return result;
}
//////////////////////////////////////////////////////////////////////////////
std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
MeshFittingAlgorithmFR::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_stiffness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFR::setup()
{
  CFAUTOTRACE;

  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;

  DataProcessingCom::setup();
  //resize and initialize the storage of the nodal stiffness
  DataHandle< CFreal > stiffness = socket_stiffness.getDataHandle();
  stiffness.resize(socket_nodes.getDataHandle().size());
  stiffness = 0.;

    // AL: this might be useless ... (@see MeshRigidMove/StdSetup.cxx)
  SubSystemStatusStack::getActive()->setMovingMesh(true);
   
  //m_lss = getMethodData().getLinearSystemSolver()[0];
  const std::string name = getMethodData().getNamespace();
  m_lss = getMethodData().getCollaborator<LinearSystemSolver>(name);  

  CFLog(VERBOSE, "MeshFittingAlgorithmFR::setup() -----=> LSS is " << m_lss->getName() << "\n");
  
  Common::SafePtr<Framework::SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();

  Common::SafePtr<FluxReconstructionSolver> frsolver = spaceMethod.d_castTo<FluxReconstructionSolver>(); //# Change here
  cf_assert(frsolver.isNotNull());
  m_frData = frsolver->getData();
  

  vector< FluxReconstructionElementData* >& frLocalData = m_frData->getFRLocalData();

 
  getMethodData().getUpdateVarSet()->setup();


  if (m_monitorVarID == std::numeric_limits<CFuint>::max() || 
    m_monitorVarID > PhysicalModelStack::getActive()->getNbEq()) {
    CFLog(WARN, "MeshFittingAlgorithmFR::setup() => monitorVarID not specified or invalid: will be set to 0\n");
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




  createGeneralConnectivityFR();    // New


  createNodalConnectivity();

  findBoundaryNodes();
  //computeMovingInBoundaryNodeNormals();   // New

}

//////////////////////////////////////////////////////////////////////////////


void MeshFittingAlgorithmFR::unsetup()
{
  CFAUTOTRACE;
  Framework::DataProcessingCom::unsetup();
}
  

//////////////////////////////////////////////////////////////////////////////


 void MeshFittingAlgorithmFR::createNodalConnectivity()
{ 
  /////////////////////////
  //2D quadrilateral mesh//
  /////////////////////////
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  // Connectivity information 2D quadrilateral
    CFuint nbPairsNodeNode = 200000; 
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeNode(nbPairsNodeNode);
  
  std::multimap<CFuint, CFuint>  mapNodeNode;
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  CFuint coupleDone [nbPairsNodeNode][2];
  for (CFuint k= 0 ; k<nbPairsNodeNode ; ++k){
    coupleDone[k][0]=0;
    coupleDone[k][1]=0;
  }
  CFuint i = 0;


  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();

    std::vector< Framework::Node*  >* m_cellNodes = currCell->getNodes();

    const CFuint nbNodes = m_cellNodes->size(); 



    const std::vector<Framework::GeometricEntity*  >* facesInCell = currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell->size(); 



      for (CFuint iFace=0; iFace<nbFaces; ++iFace){
	std::vector<Framework::Node*>* faceNodes = (*facesInCell)[iFace]->getNodes();
	CFuint nodeIDinF1 = (*faceNodes)[0]->getLocalID();
	CFuint nodeIDinF2 = (*faceNodes)[1]->getLocalID();
	  bool done1 = false;
	  bool done2 = false;

	  for (CFuint j=0; j<nbPairsNodeNode ; ++j){
	    if(coupleDone[j][0] ==nodeIDinF1 &&  coupleDone[j][1] ==nodeIDinF2){
	      done1 = true;
	    }
	    if(coupleDone[j][0] ==nodeIDinF2 &&  coupleDone[j][1] ==nodeIDinF1){
	      done2 = true;
	    }
	  }
	  if (done1==false){
	    m_mapNodeNode.insert(nodeIDinF1,nodeIDinF2);
	    coupleDone[i][0]= nodeIDinF1;
	    coupleDone[i][1]= nodeIDinF2;
	    i+=1;
	  }
	  if (done2==false){
	    m_mapNodeNode.insert(nodeIDinF2,nodeIDinF1);
	    coupleDone[i][0]= nodeIDinF2;
	    coupleDone[i][1]= nodeIDinF1;
	    i+=1;
	  }
      }
      m_mapNodeNode.sortKeys();
      
      m_cellBuilder->releaseGE();
  }
  //   }
  m_mapNodeNode1=m_mapNodeNode;

  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR();

}
//////////////////////////////////////////////////////////////////////////////*/




//////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFR::execute()  /// needs to be changed 
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithmFR::execute() => start \n");

  resizeSystemSolverToNodalData();

  computeNodeStates();

  computeSpringTruncationData(); 

  solveLinearSystem();

  // mesh node repositionning
  updateNodePositions();

  resizeSystemSolverToStateData();

  triggerRecomputeMeshData();
  CFLog(INFO, "MeshFittingAlgorithm::execute() => end \n");
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFR::computeSpringTruncationData() 
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
    //compute store the spring constants in a temporary buffer
    for(CFuint iBuffer=0; iBuffer<sizeSendBuffer && it != m_edgeGraph.end(); ++it){
      Framework::Node* node1 = (*it).firstNode;
      Framework::Node* node2 = (*it).secondNode;
      if(!isNodeLocked(node1) && node1->isParUpdatable()){
        const CFreal springConstant = computeSpringConstant(node1, node2);
        if(springConstant > 1.e-4){
          sendBuffer.push_back(springConstant);
          ++iBuffer;
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




//////////////////////////////////////////////////////////////////////////////

  CFreal MeshFittingAlgorithmFR::computeSpringConstant(const Framework::Node* const firstNode, 
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

CFreal MeshFittingAlgorithmFR::truncateSpringConstant(const CFreal springConstant){
  const CFreal maxLimit = m_springTruncationData.maxLimit;
  const CFreal minLimit = m_springTruncationData.minLimit;
  const CFreal mean     = m_springTruncationData.mean;

  const CFreal truncatedSpringConstant = std::max(std::min(maxLimit/mean, springConstant/mean), minLimit/mean );
  return truncatedSpringConstant;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFR::solveLinearSystem(){
  assembleLinearSystem();
  m_lss->solveSys();
	
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFR::resizeSystemSolverToNodalData(){
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

void MeshFittingAlgorithmFR::resizeSystemSolverToStateData(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  if (rhs.size()/totalNbEqs != states.size()) rhs.resize(states.size()*totalNbEqs); 
}



//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFR::computeNodeStates()
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
      
//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFR::assembleLockedNode(const Framework::Node* node){
     Framework::DataHandle< CFreal > stiffness = socket_stiffness.getDataHandle();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;

  //cout<< "LN: Global node ID   " << globalID << "   Local ID    " << node->getLocalID()<<endl;
   stiffness[node->getLocalID()]=100.;
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, 1.);
  }
  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = (*node)[XX+iDim];
  }
}

//////////////////////////////////////////////////////////////////////////////
 
 void MeshFittingAlgorithmFR::assembleInnerNode(const Framework::Node* node){
   SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
     getTrs("InnerCells");
   const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
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
   if( nbElemTypes==4 && nbDims==2){
     CFreal sumOffDiagonalValues = 0.;
     typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
     typedef MapNodeNode::MapIterator mapIt;
     bool found = false;
     std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
     cf_assert(found);    
      for(itN=neighboringNodes.begin(); itN != neighboringNodes.end(); ++itN){ 
       const  Framework::Node* neighborNode = *itN;
       for (mapIt it = ite.first; it != ite.second; ++it) {
            if (neighborNode->getLocalID() == it->second){
	      const CFreal springConstant =computeSpringConstant(node,neighborNode);
              const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
	      //cout<<" normalized spring constant  "<< normalizedSpringConstant << endl;
              sumOffDiagonalValues +=normalizedSpringConstant ;
              const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
	      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
  //cout<< "IN: Global node ID RG  " << rowGlobalID << "   Local ID    " << node->getLocalID()<<endl;
  //cout<< "IN: Global node ID CG  " << colGlobalID << "   Local ID    " << neighborNode->getLocalID()<<endl;
  //cout<< "idxMapping size        "<< idxMapping.getRowID(node->getLocalID()) << endl;
  //cout<< "IN: Global node ID CG  " << colGlobalID << "   Local ID    " << neighborNodes[iNode]->getLocalID()<<endl;
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
void MeshFittingAlgorithmFR::updateNodePositions () {
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::updateNodePositions()\n");  
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle < CFreal > rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
	//if (nodes[iNode]->isParUpdatable()){ 
	  
    if (nodes[iNode]->isParUpdatable()) {
      Framework::Node& currNode = *nodes[iNode];
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
        currNode[XX+iDim] = currNode[XX+iDim]*(1.-m_meshAcceleration) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration; 
      }
}
    //}
  }
  //synchronize Nodes
  nodes.beginSync();
  nodes.endSync();
}
/////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFR::assembleLinearSystem(){
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::assembleLinearSystem()\n");
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();
 
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode){
	if (!nodes[iNode]->isParUpdatable()){ 
	  //do nothing
	}
    else{
    if(isNodeMovingInBoundary(nodes[iNode])){

      assembleMovingInBoundaryNode(nodes[iNode]);
    }
    else{
      if(isNodeLocked(nodes[iNode])){
	assembleLockedNode(nodes[iNode]);
	
      }
      else{
	assembleInnerNode(nodes[iNode]);
      }
    }
  }
}
}
/////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFR::isNodeLocked( Framework::Node* node)
{
  const bool isBoundary = m_boundaryNodes.find(node) != m_boundaryNodes.end();
  return ((isBoundary)) ;
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFR::findBoundaryNodes()
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::createConnectivity()" << "\n");
  
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  std::vector<CFuint> alreadyComputed;

  for (CFuint iCell=0; iCell<nbCells; ++iCell) {
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    const std::vector<Framework::GeometricEntity*>* facesInCell = currCell->getNeighborGeos();
    Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryCell = 
      m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();

    const CFuint nbFaces = facesInCell->size();
    for (CFuint iFace=0; iFace<nbFaces; ++iFace) {
      if((*m_isFaceOnBoundaryCell)[iFace]) {
	std::vector<Framework::Node*>* faceNodes = (*facesInCell)[iFace]->getNodes();
	const CFuint nbFaceNodes = faceNodes->size();
	bool done = false;
	for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {
	  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode((*faceNodes)[iNode]);

	  for(CFuint jxx=0; jxx<alreadyComputed.size(); ++jxx){
	    if(alreadyComputed[jxx] ==  (*faceNodes)[iNode]->getLocalID()){
	      done = true;
	    }
	  }

	  if(done==false &&  neighboringNodes.size()<6 ){

	    	alreadyComputed.push_back((*faceNodes)[iNode]->getLocalID());
	    	m_boundaryNodes.insert((*faceNodes)[iNode]);
		
	  }
	}
      }
    }
    m_cellBuilder->releaseGE();
  }
}


//////////////////////////////////////////////////////////////////////////////
   
void MeshFittingAlgorithmFR::computeMovingInBoundaryNodeNormals()
{

// FB : NOT Active YET
  CFAUTOTRACE;
  const CFuint nbDim = Framework::PhysicalModelStack::getActive()->getDim();
  Framework::DataHandle<CFreal> normals = socket_normals.getDataHandle(); 
  std::vector<std::set<CFuint> > mapNodeID2BoundaryFaceIDs = getMapMovingInBoundaryNodeID2FaceID();

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
          /*if(itFaceID1 != itFaceID2){

            for(CFuint iDim=0; iDim<nbDim; ++iDim){
	      cout<< "  (*itFaceID1)                "<<(*itFaceID1)<<endl;
	      cout<< "  (*itFaceID1)*nbdims +idim   "<<(*itFaceID1)*nbDim+iDim<<endl;
	      cout<< normals.size() << endl;
              normal1[iDim] = normals[(*itFaceID1)*nbDim+iDim];
              normal2[iDim] = normals[(*itFaceID2)*nbDim+iDim];
            }
            normal1.normalize(); normal2.normalize();
            const CFreal cosAngle = MathTools::MathFunctions::innerProd(normal1, normal2);
	    cout<< " cosine angle    " << cosAngle<< endl;
            areNormalsCollinear &= (1.- std::abs(cosAngle)< 1e-5);
          }*/
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
            averageNormal[iDim] += normals[normalID*nbDim+iDim];

          }
        }
        averageNormal.normalize();
        CFLog(NOTICE,"averageNormal  "<<averageNormal<< "/n");
        m_mapNodeIDNormal.insert(addPair);
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
 
std::vector<std::set<CFuint> > MeshFittingAlgorithmFR::getMapMovingInBoundaryNodeID2FaceID()
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::computeNodes2BoundaryNormals()" << "\n");
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();
  std::vector<std::set<CFuint> > mapBoundaryNodeID2FaceIDs(nbNodes);
  for(CFuint iTrs=0; iTrs<m_unlockedBoundaryTRSs.size();++iTrs){
    FaceToCellGEBuilder::GeoData& facesData = m_faceBuilder->getDataGE();
    m_faceBuilder->getDataGE().isBoundary = true;
    Common::SafePtr<Framework::TopologicalRegionSet> wallFaces =
    Framework::MeshDataStack::getActive()->getTrs( m_unlockedBoundaryTRSs[iTrs] );
    facesData.facesTRS = wallFaces;
    const CFuint nbFaces = wallFaces->getLocalNbGeoEnts();

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
          (mapBoundaryNodeID2FaceIDs[nodeID]).insert(faceID);
        }
      }
      m_faceBuilder->releaseGE();
    }
  }
  return mapBoundaryNodeID2FaceIDs;
}

//////////////////////////////////////////////////////////////////////////////

bool MeshFittingAlgorithmFR::isNodeMovingInBoundary(Framework::Node* node){
  
  return (m_mapNodeIDNormal.find(node->getLocalID()) != m_mapNodeIDNormal.end());
}


//////////////////////////////////////////////////////////////////////////////

 void MeshFittingAlgorithmFR::createGeneralConnectivityFR()
{ 
  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR(); 
  }


//////////////////////////////////////////////////////////////////////////////
 
void MeshFittingAlgorithmFR::triggerRecomputeMeshData() {
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
