#include "Common/PE.hh"
#include "Common/EventHandler.hh"
#include "MathTools/LeastSquaresSolver.hh"
#include "MathTools/MathFunctions.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BaseTerm.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <iostream>
#include <limits>

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshFittingAlgorithm, 
		      DataProcessingData, 
		      FiniteVolumeModule>
meshFittingAlgorithmProvider("MeshFittingAlgorithm");

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::defineConfigOptions(Config::OptionList& options)
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

//////////////////////////////////////////////////////////////////////////////

MeshFittingAlgorithm::MeshFittingAlgorithm(const std::string& name) :
  Framework::DataProcessingCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_gstates("gstates"),
  socket_normals("normals"),
  socket_rhs("rhs"), 
  socket_wallDistance("wallDistance",false),
  m_wallDistance(CFNULL),
  m_lss(CFNULL),
  m_fvmccData(CFNULL),
  m_geoBuilder(),
  m_faceTRSBuilder(),
  m_pdata()
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

MeshFittingAlgorithm::~MeshFittingAlgorithm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::configure(Config::ConfigArgs& args)
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
MeshFittingAlgorithm::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);
  result.push_back(&socket_rhs);
  result.push_back(&socket_wallDistance);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
MeshFittingAlgorithm::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::setup()
{
  CFAUTOTRACE;

  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  
  // AL: this might be useless ... (@see MeshRigidMove/StdSetup.cxx)
  SubSystemStatusStack::getActive()->setMovingMesh(true);
  
  // get the linear system associated to this object and set it up
  const std::string name = getMethodData().getNamespace();
  m_lss = getMethodData().getCollaborator<LinearSystemSolver>(name);  
  
  CFLog(VERBOSE, "MeshFittingAlgorithm::setup() => LSS is " << m_lss->getName() << "\n");
  
  Common::SafePtr<Framework::SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  Common::SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());
  m_fvmccData = fvmcc->getData();
  
  /// setup geobuilder
  m_geoBuilder.setup();
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  m_faceTRSBuilder.setup();
  m_faceTRSBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);

  getMethodData().getUpdateVarSet()->setup();
  
  // resize the physical data array
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(m_pdata);
  
  m_state.reset(new State());
  
  if (m_monitorVarID == std::numeric_limits<CFuint>::max() || 
      m_monitorVarID > PhysicalModelStack::getActive()->getNbEq()) {
    CFLog(WARN, "MeshFittingAlgorithm::setup() => monitorVarID not specified or invalid: will be set to 0\n");
    m_monitorVarID = 0;
  }
  
  if (m_monitorPhysVarID != std::numeric_limits<CFuint>::max() && 
      m_monitorPhysVarID >= m_pdata.size()) {
    CFLog(WARN, "MeshFittingAlgorithm::setup() => monitorPhysVarID not valid: monitorVarID will be used instead\n");
    // m_monitorPhysVarID is reset to default
    m_monitorPhysVarID = std::numeric_limits<CFuint>::max();
  }
  
  const std::string namespaceName = MeshDataStack::getActive()->getPrimaryNamespace();
  const string wallDistanceDataHandleName = namespaceName + "_wallDistance";
  const bool wallDistanceExists = MeshDataStack::getActive()->
    getDataStorage()->checkData(wallDistanceDataHandleName);
  if (wallDistanceExists) {
    m_wallDistance = MeshDataStack::getActive()->getDataStorage()->
      getData<CFreal>(wallDistanceDataHandleName);
    cf_assert(m_wallDistance.getLocalArray() != CFNULL);
    cf_assert(m_wallDistance.size() > 0);
  }
  
  createNodalConnectivity();
  findBoundaryNodes();
  computeMovingInBoundaryNodeNormals();
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::createNodalConnectivity()
{
  m_edgeGraph.setNodeDataSocket(socket_nodes);
  m_edgeGraph.computeConnectivity(); 
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::findBoundaryNodes()
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::createConnectivity()" << "\n");
  
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  
  for (CFuint iCell=0; iCell<nbCells; ++iCell) {
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size();
    
    for (CFuint iFace=0; iFace<nbFaces; ++iFace) {
      std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
      const CFuint nbFaceNodes = faceNodes.size();

      if (facesInCell[iFace]->getState(1)->isGhost()) {
        for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {
          m_boundaryNodes.insert(faceNodes[iNode]);
        }
      }
    }
    m_geoBuilder.releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////
   
void MeshFittingAlgorithm::computeMovingInBoundaryNodeNormals()
{
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
          if(itFaceID1 != itFaceID2){
            for(CFuint iDim=0; iDim<nbDim; ++iDim){
              normal1[iDim] = normals[(*itFaceID1)*nbDim+iDim];
              normal2[iDim] = normals[(*itFaceID2)*nbDim+iDim];
            }
            normal1.normalize(); normal2.normalize();
            const CFreal cosAngle = MathTools::MathFunctions::innerProd(normal1, normal2);
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
            averageNormal[iDim] += normals[normalID*nbDim+iDim];
          }
        }
        averageNormal.normalize();
        m_mapNodeIDNormal.insert(addPair);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
 
std::vector<std::set<CFuint> > MeshFittingAlgorithm::getMapMovingInBoundaryNodeID2FaceID()
{
  CFAUTOTRACE;
  CFLogDebugMin("MeshFittingAlgorithm::computeNodes2BoundaryNormals()" << "\n");
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  const CFuint nbNodes = nodes.size();
  std::vector<std::set<CFuint> > mapBoundaryNodeID2FaceIDs(nbNodes);
  for(CFuint iTrs=0; iTrs<m_unlockedBoundaryTRSs.size();++iTrs){
    Framework::FaceTrsGeoBuilder::GeoData& facesData = m_faceTRSBuilder.getDataGE();
    m_faceTRSBuilder.getDataGE().isBFace = true;
    Common::SafePtr<Framework::TopologicalRegionSet> wallFaces =
    Framework::MeshDataStack::getActive()->getTrs( m_unlockedBoundaryTRSs[iTrs] );
    facesData.trs = wallFaces;
    const CFuint nbFaces = wallFaces->getLocalNbGeoEnts();
    for(CFuint i=0; i<nbFaces; ++i){
      facesData.idx = i;
      const Framework::GeometricEntity *const face = m_faceTRSBuilder.buildGE();
      const std::vector<Framework::Node*> faceNodes = face->getNodes();
      const CFuint faceID = face->getID();
      const CFuint nbFaceNodes = faceNodes.size();
      if (face->getState(1)->isGhost()) {
        for(CFuint iNode=0; iNode<nbFaceNodes; ++iNode) {
          const CFuint nodeID = faceNodes[iNode]->getLocalID();
          (mapBoundaryNodeID2FaceIDs[nodeID]).insert(faceID);
        }
      }
      m_faceTRSBuilder.releaseGE();
    }
  }
  return mapBoundaryNodeID2FaceIDs;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::execute()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::execute() => start \n");
  
  resizeSystemSolverToNodalData();
  computeSpringTruncationData();
  solveLinearSystem();
  updateNodePositions();
  resizeSystemSolverToStateData();
  triggerRecomputeMeshData();
  
  CFLog(INFO, "MeshFittingAlgorithm::execute() => end \n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::computeSpringTruncationData() 
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
  
  SimpleEdgeGraph::iterator it = m_edgeGraph.begin();
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
  
CFreal MeshFittingAlgorithm::computeSpringConstant(const Framework::Node* const firstNode, 
                                                          const Framework::Node* const secondNode) 
{
  CFAUTOTRACE;
  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
    const CFreal firstNodeValue  = nodalStates[firstNode->getLocalID()] [m_monitorVarID];
    const CFreal secondNodeValue = nodalStates[secondNode->getLocalID()][m_monitorVarID];
    return std::abs(secondNodeValue - firstNodeValue);
  }
  
  cf_assert(m_monitorPhysVarID < m_pdata.size());
  // physical data arrays are computed on-the-fly from given nodal states 
  m_state->copyData(nodalStates[firstNode->getLocalID()]);
  getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
  const CFreal firstNodeValue  = m_pdata[m_monitorPhysVarID];
  m_state->copyData(nodalStates[secondNode->getLocalID()]);
  getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
  const CFreal secondNodeValue  = m_pdata[m_monitorPhysVarID];
  return std::abs(secondNodeValue - firstNodeValue);
}

//////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithm::truncateSpringConstant(const CFreal springConstant){
  const CFreal maxLimit = m_springTruncationData.maxLimit;
  const CFreal minLimit = m_springTruncationData.minLimit;
  const CFreal mean     = m_springTruncationData.mean;

  const CFreal truncatedSpringConstant = std::max(std::min(maxLimit/mean, springConstant/mean), minLimit/mean );
  return truncatedSpringConstant;
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::solveLinearSystem(){
  assembleLinearSystem();
  m_lss->solveSys();
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::resizeSystemSolverToNodalData(){
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint rhsSize = std::max(nodes.size()*totalNbEqs,states.size()*totalNbEqs);
  if (rhs.size()!=rhsSize) rhs.resize(rhsSize);
  rhs = 0.;
  jacobMatrix->resetToZeroEntries();
} 
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::resizeSystemSolverToStateData(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  if (rhs.size()/totalNbEqs != states.size()) rhs.resize(states.size()*totalNbEqs); 
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::assembleLinearSystem(){
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::assembleLinearSystem()\n");

  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  //const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim(); 
  const CFuint nbNodes = nodes.size();
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
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

//////////////////////////////////////////////////////////////////////////////

bool MeshFittingAlgorithm::isNodeMovingInBoundary(Framework::Node* node){
  return (m_mapNodeIDNormal.find(node->getLocalID()) != m_mapNodeIDNormal.end());
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::assembleMovingInBoundaryNode(const Framework::Node* node){
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  const RealVector& nodeNormal = (m_mapNodeIDNormal[node->getLocalID()]);

  //Choose the dependent dimension and the free dimension 
  //by doing so we avoid problems with the pivoting of the LU decomposition
  
  //Dependent dim is the dimension with highest normal value
  CFuint dependentDim = 0;
  for (CFuint i=1; i<nbDims; ++i){
    dependentDim = (std::abs(nodeNormal[i]) > std::abs(nodeNormal[dependentDim])) ? i : dependentDim ;
  }
  //Free dims are the others 
  std::vector<CFuint> freeDims;
  for (CFuint iDim=0; iDim<nbDims; ++iDim){
    if (iDim!=dependentDim){
       freeDims.push_back(iDim);
    }
  }

  //Dependent dimension, following the line equation a*y + b*x + c*y= d,
  // where y is the dep dim, x is the free dim, a and b are the normals and d is the y-intersect
  const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
  CFreal y_intersect = nodeNormal[dependentDim]*(*node)[XX+dependentDim];
  for(CFuint iFreeDim=0; iFreeDim<freeDims.size();++iFreeDim){
    y_intersect += nodeNormal[freeDims[iFreeDim]]*(*node)[XX+freeDims[iFreeDim]];
  }
  jacobMatrix->addValue(globalID+dependentDim, globalID+dependentDim, nodeNormal[dependentDim]);
  for(CFuint iFreeDim=0; iFreeDim<freeDims.size();++iFreeDim){
    jacobMatrix->addValue(globalID+dependentDim, globalID+freeDims[iFreeDim], nodeNormal[freeDims[iFreeDim]]);
  }
  rhs[node->getLocalID()*totalNbEqs+XX+dependentDim] = y_intersect;

  //Free dimensions, acting like an inner Node with projected spring constants
  CFreal sumOffDiagonalValues = 0.;
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
  std::vector<Framework::Node*>::const_iterator it;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    const Framework::Node* neighborNode = *it;
    const CFreal springConstant = computeSpringConstant(node, neighborNode);

    RealVector springDirection(nbDims);
    for (CFuint iDim=0; iDim<nbDims; ++iDim){
      springDirection[iDim] = (*neighborNode)[XX+iDim] - (*node)[XX+iDim] ;
    }
    springDirection.normalize();
    RealVector projectedSpringDirection = springDirection - 
      MathTools::MathFunctions::innerProd(springDirection, nodeNormal)*nodeNormal;

    const CFreal dotProduct = MathTools::MathFunctions::innerProd(projectedSpringDirection, springDirection);
    const CFreal normalizedSpringConstant = truncateSpringConstant(springConstant)*std::abs(dotProduct);
    sumOffDiagonalValues += normalizedSpringConstant;
    
    const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
    const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
    for(CFuint iFreeDim=0; iFreeDim<freeDims.size(); ++iFreeDim){
      jacobMatrix->addValue(rowGlobalID+freeDims[iFreeDim], colGlobalID+freeDims[iFreeDim], normalizedSpringConstant);
    }
  }
  for(CFuint iFreeDim=0; iFreeDim<freeDims.size(); ++iFreeDim){
    CFreal diagValue = -sumOffDiagonalValues;
    jacobMatrix->addValue(globalID+freeDims[iFreeDim], globalID+freeDims[iFreeDim], diagValue);
    rhs[node->getLocalID()*totalNbEqs+XX+freeDims[iFreeDim]] = m_equilibriumSpringLength*sumOffDiagonalValues;
  }
}

//////////////////////////////////////////////////////////////////////////////

bool MeshFittingAlgorithm::isNodeLocked(Framework::Node* node)
{
  const bool isBoundary = m_boundaryNodes.find(node) != m_boundaryNodes.end();
  return (isBoundary);
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::assembleLockedNode(const Framework::Node* node){

  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, 1.);
  }

  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = (*node)[XX+iDim];
  }
}

//////////////////////////////////////////////////////////////////////////////
 
void MeshFittingAlgorithm::assembleInnerNode(const Framework::Node* node){
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  CFreal sumOffDiagonalValues = 0.;
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
  std::vector<Framework::Node*>::const_iterator it;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    const Framework::Node* neighborNode = *it;
    const CFreal springConstant = computeSpringConstant(node, neighborNode);
    const CFreal normalizedSpringConstant = truncateSpringConstant(springConstant);
    sumOffDiagonalValues += normalizedSpringConstant;
    const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
    const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
    for(CFuint iDim=0; iDim<nbDims; ++iDim){
      jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim, normalizedSpringConstant);
    }
  }
  //Diagonal value 
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue = -sumOffDiagonalValues;
    const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
  }

  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*sumOffDiagonalValues;
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::updateNodePositions() {
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::updateNodePositions()\n");
  
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle < CFreal > rhs = socket_rhs.getDataHandle();
  
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {
    if (nodes[iNode]->isParUpdatable()) {
      Framework::Node& currNode = *nodes[iNode];
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
        currNode[XX+iDim] = currNode[XX+iDim]*(1.-m_meshAcceleration) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration;
      }
    }
  }
  //synchronize Nodes
  nodes.beginSync();
  nodes.endSync();
}

//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::triggerRecomputeMeshData() {
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
      
//////////////////////////////////////////////////////////////////////////////


    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
