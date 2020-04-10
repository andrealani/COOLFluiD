#include "Common/PE.hh"
#include "Common/EventHandler.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"
//#include "MathTools/LeastSquaresSolver.hh"



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
MeshFittingAlgorithmFRQ2FluxReconstructionProvider("MeshFittingAlgorithmFRQ2");

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
}

///////////////////////////////////////////////////////////////////

MeshFittingAlgorithmFRQ2::MeshFittingAlgorithmFRQ2(const std::string& name) :

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
  //m_nbOfNeighborCellsToaNode(CFNULL)
  //m_boundaryNodes(CFNULL),
  //m_CenterNodes(CFNULL)
  
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

  //resize and initialize the storage of the nodal stiffness
  DataHandle< CFreal > stiffness = socket_stiffness.getDataHandle();
  stiffness.resize(socket_nodes.getDataHandle().size());
  stiffness = 0.;

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();

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

  m_edgeGraph.setNodeDataSocketFR(socket_nodes);
  m_edgeGraph.computeConnectivityFR();

  std::vector<CFuint >  nbOfNeighborCellsToaNode = m_edgeGraph.getNbOfNeighborCells();
  m_nbOfNeighborCellsToaNode = nbOfNeighborCellsToaNode;

  createGeneralConnectivityFR();    

  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();

  findBoundaryNodes();

  findCenterNodes();


  nbOfConnectedFacesToaNode();

  createNodalConnectivity();
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle(); 
  CFuint nbNodes = nodes.size();
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();



  // get the coefs for extrapolation of the states to the flx pnts
  //m_solPolyValsAtNodes = frLocalData[0]->getCoefSolPolyInNodes(); // OLD only 4 corners

  /*Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  CFuint i=0;

  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    
    std::vector< Framework::Node*  >* m_cellNodes = currCell->getNodes();
    
    const CFuint nbNodes = m_cellNodes->size(); 
    for(CFuint iNode=0; iNode<nbNodes; ++iNode){
      Framework::Node * currNode = (*m_cellNodes)[iNode];


      i=i+1;
    }
    m_cellBuilder->releaseGE();
  }*/
    std::vector< RealVector > vec_cor_ALL; 

  for (CFuint iNode=0; iNode<nbNodes; ++iNode){
    Framework::Node * currNode = nodes[iNode];
    RealVector vec_cor; vec_cor.resize(nbDims);

    for(CFuint iDim=0; iDim<nbDims; ++iDim ){
      vec_cor[iDim] = (*currNode)[XX+iDim];
    }
    vec_cor_ALL.push_back(vec_cor);
  }
    m_solPolyValsAtNodes = frLocalData[0]->getNodePolyValsAtPnt(vec_cor_ALL);
    m_nbrNodesElem = m_solPolyValsAtNodes.size();


}

//////////////////////////////////////////////////////////////////////////////


void MeshFittingAlgorithmFRQ2::unsetup()
{
  CFAUTOTRACE;
  Framework::DataProcessingCom::unsetup();
}
  

//////////////////////////////////////////////////////////////////////////////


 void MeshFittingAlgorithmFRQ2::createNodalConnectivity()
{ 
  /////////////////////////
  //2D quadrilateral mesh//
  /////////////////////////

  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle(); 

  // Connectivity information 2D quadrilateral
  CFuint nbPairsNodeNode = 20000; 
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  Common::CFMultiMap<CFuint, CFuint>  m_mapNodeNode(nbPairsNodeNode);
  
  std::multimap<CFuint, CFuint>  mapNodeNode;
  std::vector< bool > nodeDone; nodeDone.resize(nodes.size());
  std::vector< CFuint > BCnodeDone; BCnodeDone.resize(nodes.size());  // ount conectivity to each node 


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
            for(CFuint iNode=0; iNode< nbNodes; ++ iNode){
              if (isInSideCell((*m_cellNodes)[iNode]) && alreadyComputed == false){
                m_mapNodeNode.insert(currNode.getLocalID(),(*m_cellNodes)[iNode]->getLocalID()); // Adding the Center point 
                BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
                alreadyComputed = true;
              }
            }
          }
          nodeDone[currNode.getLocalID()] = true;
        }




        if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] == 2 && isBoundaryNode(&currNode)){ // Our node: Boundary intersection Point
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
            BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
          }
        }

      if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] == 1 && isBoundaryNode(&currNode) && nbOfConnectedFaces[currNode.getLocalID()] == 2){ // Our node: Corner Point
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
        }



        if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] == 2 && isBoundaryNode(&currNode)==false){ // Our node: Middle point between 2 cells            
            for(CFuint iNode=0; iNode< nbNodes; ++ iNode){
              if (isInSideCell((*m_cellNodes)[iNode])){
                m_mapNodeNode.insert(currNode.getLocalID(),(*m_cellNodes)[iNode]->getLocalID()); // Adding the Center point 
                BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;              
              }
            }
            for (CFuint iNode=0; iNode<nbNodesinF; ++iNode){
              if((faceNodes)[iNode]->getLocalID() != currNode.getLocalID()   && coupleDone[currNode.getLocalID()][(faceNodes)[iNode]->getLocalID()]==false){
                m_mapNodeNode.insert(currNode.getLocalID(), (faceNodes)[iNode]->getLocalID()); // Adding the 2 end points 
                coupleDone[currNode.getLocalID()][(faceNodes)[iNode]->getLocalID()]=true;  // Cheking if they are done
                BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
              }
            }
          nodeDone[currNode.getLocalID()] = true;
        }




        if (m_nbOfNeighborCellsToaNode[currNode.getLocalID()] == 4){ // Our node: intersection between 2 edges
          
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
          BCnodeDone[currNode.getLocalID()]=BCnodeDone[currNode.getLocalID()]+1;
          nodeDone[currNode.getLocalID()] = true;
          coupleDone[currNode.getLocalID()][closestNode]=true;

          }
        }
      }
    }
    m_cellBuilder->releaseGE();
  }

  for (CFuint iNode=0; iNode<nodes.size(); ++iNode){
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
  }




  m_mapNodeNode.sortKeys();
  //   }
  m_mapNodeNode1=m_mapNodeNode;
          //cout<<" here10 "<<endl;
          //cout<<" Finish Connectvity "<<endl;
}

//////////////////////////////////////////////////////////////////////////////*/




//////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::execute()  /// needs to be changed 
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithmFRQ2::execute() => start \n");

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
  
  SimpleEdgeGraphFRQ2::iterator it = m_edgeGraph.begin();
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
  ////cout<<" -----------------state size " << states.size()<< endl;
  ////cout<<" -----------------RHS size " << rhsSize<< endl;
  ////cout<<" -----------------totalNbEqs size " << totalNbEqs<< endl;
  ////cout<<" -----------------nodes size " << nodes.size()<< endl;

} 
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::resizeSystemSolverToStateData(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  if (rhs.size()/totalNbEqs != states.size()) rhs.resize(states.size()*totalNbEqs); 
}



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

////cout<<" nodal states size()  "<< nodalStates[10].size() << endl;
////cout<<"  m_solPolyValsAtNodes "<<m_solPolyValsAtNodes[10].size()<<endl;
  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();
    
    std::vector< Framework::Node*  >* m_cellNodes = currCell->getNodes();
    const CFuint nbNodes = m_cellNodes->size(); 
    ////cout<<"   nbNodes   "<< nbNodes << endl;
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
        ////cout<< "iNode  "<< iNode << "   iSol  " << iSol << "   iVar  "<< iVar<<endl;
        ////cout<<" (m_solPolyValsAtNodes)[nodes[iNode]->getLocalID()][iSol]   "  << (m_solPolyValsAtNodes)[(*m_cellNodes)[iNode]->getLocalID()][iSol] << endl;
        ////cout<<" (*(*m_cellStates)[iSol])[iVar]   "  << (*(*m_cellStates)[iSol])[iVar] << endl;
	      nodalStates[(*m_cellNodes)[iNode]->getLocalID()][iVar] += (m_solPolyValsAtNodes)[(*m_cellNodes)[iNode]->getLocalID()][iSol]*(*(*m_cellStates)[iSol])[iVar];
	      //nodalStates[(*m_cellNodes)[iNode]->getLocalID()][iVar] += (*m_solPolyValsAtNodes)[iNode][iSol]*(*(*m_cellStates)[iSol])[iVar];
        ////cout<<" nodalStates   "  << nodalStates[(*m_cellNodes)[iNode]->getLocalID()][iVar] << endl;
        ////cout<<"  "<<endl;

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
    //cout<<" counterICI " << counter[(nodes)[iNode]->getLocalID()]<< endl;
  }
}*/
  for (CFuint iNode=0 ; iNode<nodalStates.size() ; ++iNode){


    nodalStates[(nodes)[iNode]->getLocalID()] = nodalStates[(nodes)[iNode]->getLocalID()]/(m_nbOfNeighborCellsToaNode[(nodes)[iNode]->getLocalID()]*m_order);//counter[(nodes)[iNode]->getLocalID()];
    //cout << "x dim" << (*nodes[iNode])[0] << "ydim" << (*nodes[iNode])[1] << endl;
    //cout<<" nodalStates   "  << nodalStates[(nodes)[iNode]->getLocalID()][0] <<"   " << nodalStates[(nodes)[iNode]->getLocalID()][1]  << "    " << nodalStates[(nodes)[iNode]->getLocalID()][2] <<"   " << nodalStates[(nodes)[iNode]->getLocalID()][3]<<endl;
    //cout<<" counter " << m_nbOfNeighborCellsToaNode[(nodes)[iNode]->getLocalID()]<< endl;
    //cout<<"   " << endl;    
////cout<<" nodalStates   "  << nodalStates[i][0] <<"   " << nodalStates[i][1]  << "    " << nodalStates[i][2] <<"   " << nodalStates[i][3]<<endl;
    ////cout<<" nodalStates   "  << nodalStates[i][0] <<"   " << nodalStates[i][1]  << "    " << nodalStates[i][2] <<"   " << nodalStates[i][3]<<endl;
  }
}


//////////////////////////////////////////////////////////////////////////////
      
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
        const  Framework::Node* neighborNode = *itN;
        for (mapIt it = ite.first; it != ite.second; ++it) {
          if(neighborNode->getLocalID() == it->second){
	          const CFreal springConstant =computeSpringConstant(node,neighborNode);
            const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
	          // //cout<<" normalized spring constant  "<< normalizedSpringConstant << endl;
            sumOffDiagonalValues +=normalizedSpringConstant ;
            const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
	          const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
            ////cout<< "IN: Global node ID RG  " << rowGlobalID << "   Local ID    " << node->getLocalID()<<endl;
            ////cout<< "IN: Global node ID CG  " << colGlobalID << "   Local ID    " << neighborNode->getLocalID()<<endl;
            ////cout<< "idxMapping size        "<< idxMapping.getRowID(node->getLocalID()) << endl;
            ////cout<< "IN: Global node ID CG  " << colGlobalID << "   Local ID    " << neighborNode->getLocalID()<<endl;
	          for(CFuint iDim=0; iDim<nbDims; ++iDim){
	            jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,normalizedSpringConstant);
            }
	          stiffness[node->getLocalID()]=normalizedSpringConstant;
          }
        }
      }
     for(CFuint iDim=0; iDim<nbDims; ++iDim){     
       const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;
       ////cout<< "diagvalues  "<<globalID << endl;
       jacobMatrix->addValue(globalID+iDim, globalID+iDim, -sumOffDiagonalValues);
     }

     //Right hand side
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
       ////cout<<" RHS Filling " << endl;
       rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues);
     }
     //jacobMatrix->printToScreen();

 }
/////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithmFRQ2::updateNodePositions () {
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
void MeshFittingAlgorithmFRQ2::assembleLinearSystem(){
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::assembleLinearSystem()\n");
  Framework::DataHandle<Framework::Node* , Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();
 
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode){
	if (!nodes[iNode]->isParUpdatable()){ 
	  //do nothing
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
//}
/////////////////////////////////////////////////////////////////////////////
bool MeshFittingAlgorithmFRQ2::isNodeLocked( Framework::Node* node)
{
  const bool isBoundary = m_boundaryNodes.find(node) != m_boundaryNodes.end();
  return ((isBoundary)) ;
}
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::findBoundaryNodes()
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::createConnectivity()" << "\n");
  
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();

  for (CFuint iCell=0; iCell<nbCells; ++iCell) {
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_cellBuilder->buildGE();

    const std::vector<Framework::GeometricEntity*  >& facesInCell = *currCell->getNeighborGeos();
    Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryCell = 
      m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();

    const CFuint nbFaces = facesInCell.size();
    for (CFuint iFace=0; iFace<nbFaces; ++iFace) {
      if((*m_isFaceOnBoundaryCell)[iFace]) {
	      std::vector<Framework::Node* >& faceNodes = *facesInCell[iFace]->getNodes();
	      const CFuint nbFaceNodes = faceNodes.size();
	      for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {

          m_boundaryNodes.insert(faceNodes[iNode]);
	      } 
	    }
    }
    m_cellBuilder->releaseGE();
  }

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
	    	    
          }
        }
      }
    }
  m_cellBuilder->releaseGE();
  }
  for (CFuint iNode=0; iNode<nodes.size(); ++iNode){
    if(m_inFace[nodes[iNode]->getLocalID()] == false){
        m_CenterNodes.insert(nodes[iNode]);
    }


  }



}


//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithmFRQ2::nbOfConnectedFacesToaNode()
{

  //FB: Only for the corner node (2 Faces) and middle boundary point (1 Face), which is the goal for this function 
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
bool MeshFittingAlgorithmFRQ2::isInSideCell(Framework::Node* node){ 
  const bool isinsideCell = m_CenterNodes.find(node) != m_CenterNodes.end();
  return (isinsideCell) ;
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
  ////cout<<" neighbor node size "<< neighboringNodes.size()<< endl;
   // }
  }


//////////////////////////////////////////////////////////////////////////////
 
void MeshFittingAlgorithmFRQ2::triggerRecomputeMeshData() {
  std::string msg;
  Common::SafePtr<Common::EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  const std::string ssname = Framework::SubSystemStatusStack::getCurrentName();   
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MESHADAPTER_AFTERMESHUPDATE"), msg );
  /*//cout << "A0" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_UNSETUP"), msg );
  ////cout << "A1" << endl;
 // event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_UNPLUGSOCKETS"), msg );
  ////cout << "A2" << endl;
  // event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_DESTROYSUBSYSTEM"), msg ); 
 // //cout << "A3" << endl;
  // event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_BUILDSUBSYSTEM"), msg );     
 ////cout << "A4" << endl;
  // event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_PLUGSOCKETS"), msg );
  //cout << "A5" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_BUILDMESHDATA"), msg );
  //cout << "A6" << endl;
  event_handler->call_signal (event_handler->key(ssname, "CF_ON_MAESTRO_SETUP"), msg );
  //cout << "A7" << endl;*/
}
    }/// Namespace FR
  /// Namespace COOLFluiD
  
}
