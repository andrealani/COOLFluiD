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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <iostream>
#include <limits>
#include <vector>
#include <typeinfo>

#include "FiniteVolume/CellData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/MeshFittingAlgorithm.hh"
#include "MeshTools/ComputeWallDistanceVector2CCMPI.hh"
#include "MeshTools/MeshToolsFVM.hh"

#include "Common/CFLog.hh"
#include "Common/CFMultiMap.hh"
#include <math.h>
#include <cmath>
#include <fstream>


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshFittingAlgorithm, 
		      DataProcessingData, 
		      FiniteVolumeModule>
meshFittingAlgorithmProvider("MeshFittingAlgorithm");
      ///////////////////////////////////////////////////////////////////////


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
  options.addConfigOption< CFreal >("tolerance","tolerance on the mesh mouvement between 2 mesh fitting process");
  // follow the MQI or the mesh stiffness :
  // 0 : deactivated
  // 1 : stiffness
  // 2 : MQI radius triangular 
  // 3 : MQI Aspect Ratio quadrilateral
  // 4 : MQI Skewness quadrilateral 
  // 5 : MQI radius sphere
  options.addConfigOption< CFuint >("MQIvalue","choosing the MQI state variable to follow");
  options.addConfigOption< CFreal >("AcceptableDistance","Distance from user-defined boundary");
  options.addConfigOption< bool   >("ThetaMid","Semi torsional Sping analogy for 2D quadrilateral mesh based on the middle egde-facing  angle (true) or the 3 edge-facing angles (false) ");
  options.addConfigOption< bool   >("InterpolateState"," State Interplation to dissociate the nodal movement and the solution in each CC. ");


  options.addConfigOption< bool >("smoothSpringNetwork","smooth the spring network");
  options.addConfigOption< bool >("smoothNodalDisp","smooth the nodal displacememt");




}

//////////////////////////////////////////////////////////////////////////////

MeshFittingAlgorithm::MeshFittingAlgorithm(const std::string& name) :
  Framework::DataProcessingCom(name),
  socket_stiffness("stiffness"),
  socket_iradius("iradius"),
  socket_skewness("skewness"),
  socket_AR("AR"),
  socket_isphere("isphere"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_gstates("gstates"),
  socket_normals("normals"),
  socket_rhs("rhs"), 
  socket_wallDistance("wallDistance",false),
  socket_nodeisAD("nodeisAD",false),
  socket_nodeDistance("nodeDistance", false),
  socket_stencil("stencil"),
  socket_relativeError("relativeError"),
  m_wallDistance(CFNULL),
  m_lss(CFNULL),
  m_fvmccData(CFNULL),
  m_geoBuilder(),
  m_faceTRSBuilder(),
  m_pdata(),
  m_mapNodeRadius1(),
  m_mapNodeNState1(),
  m_mapNodeSkew1(),
  m_mapNodeTS1(),
  m_mapNodeAR1(),
  oldStates(),
  oldCoordinates(),
  m_erreurG()

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

  m_tolerance = 0.;
  this->setParameter("tolerance", &m_tolerance);

  m_MQIvalue = 0;
  this->setParameter("MQIvalue", &m_MQIvalue);

  m_acceptableDistance = 0.;
  this->setParameter("AcceptableDistance",&m_acceptableDistance);

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
  cf_assert(m_tolerance >= 0. && m_tolerance <= 1.);
  cf_assert(m_MQIvalue >= 0);
  cf_assert(m_acceptableDistance >= 0);
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
  result.push_back(&socket_nodeisAD);  
  result.push_back(&socket_nodeDistance);
  result.push_back(&socket_stencil);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
MeshFittingAlgorithm::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_stiffness);
  result.push_back(&socket_iradius);
  result.push_back(&socket_skewness);
  result.push_back(&socket_AR);
  result.push_back(&socket_isphere);
  result.push_back(&socket_relativeError);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

 void MeshFittingAlgorithm::setup()
{
  CFAUTOTRACE;

  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  
  DataProcessingCom::setup();

  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();

  //resize and initialize the storage of the relative Eroor
  DataHandle<CFreal> relativeError = socket_relativeError.getDataHandle();
  relativeError.resize(socket_nodes.getDataHandle().size());
  relativeError = 1.;

  //resize and initialize the storage of the nodal stiffness
  DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
  stiffness.resize(socket_nodes.getDataHandle().size());
  stiffness = 0.;

  //resize and initialize the storage of the incerted radius
  DataHandle<CFreal> iradius = socket_iradius.getDataHandle();
  iradius.resize(socket_nodes.getDataHandle().size());
  iradius = 0.;

  //resize and initialize the storage of the incerted radius
  DataHandle<CFreal> skewness = socket_skewness.getDataHandle();
  skewness.resize(socket_nodes.getDataHandle().size());
  skewness = 0.;


  //resize and initialize the storage of the incerted radius
  DataHandle<CFreal> AR = socket_AR.getDataHandle();
  AR.resize(socket_nodes.getDataHandle().size());
  AR = 0.;

  //resize and initialize the storage of the incerted radius
  DataHandle<CFreal> isphere = socket_isphere.getDataHandle();
  isphere.resize(socket_nodes.getDataHandle().size());
  isphere = 0.;


  oldStates.resize(nbEqs*socket_states.getDataHandle().size());
  oldStates = 0.;
  oldCoordinates.resize(nbDims*socket_states.getDataHandle().size());
  oldCoordinates = 0.;

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
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
    getTrs("InnerCells");
  const CFuint nbElemTypes = cells->getNbNodesInGeo(0);

  ////////////////////////////////////////
  // 2D triangular mesh                 //
  // MQI based radius of inserted circle//
  ////////////////////////////////////////

  if (nbElemTypes==3 && nbDims==2 && m_MQIvalue==2){
    Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
    CFuint nbPairsNodeNode = 200000;
    typedef CFMultiMap<CFuint, CFreal> MapNodeRadius;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeRadius(nbPairsNodeNode);
    typedef CFMultiMap<CFuint, CFreal> MapNodeNState;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNState(nodes.size());
    Common::SafePtr<Framework::TopologicalRegionSet> cells =
      Framework::MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    geoData.trs = cells;
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
      geoData.idx = iCell;
      Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
      Framework::Node * thirdNode;
      const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
      const CFuint nbFaces = facesInCell.size(); 
      RealVector res(3);res=0.;
      std::vector<Framework::Node*>& faceNodes1 = *facesInCell[0]->getNodes();
      bool found=false;
      RealVector v1(3); v1=0.;
      RealVector v2(3); v2=0.;
      RealVector v3(3); v3=0.;
      
      //First order simulation v1state = user-defined
      //  variable to be putted w.r.t the monitor variable
      //  initial state  
      
      CFreal v1state=1.;
      
      ////////////////////////////////////////////////////
      //Starting from second order v1state = NodalStates//
      ////////////////////////////////////////////////////
      /*  if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
	  v1state= (*states[iCell])[m_monitorVarID];
	  }
	  else{
	  cf_assert(m_monitorPhysVarID < m_pdata.size());
	  m_state->copyData(*states[iCell]);
	  getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
	  v1state  = m_pdata[m_monitorPhysVarID];
	  }*/
      //////////////////////
      //End of second order/
      //////////////////////
      
      for (CFuint iFace=1; iFace<nbFaces; ++iFace){
	while(found==false){
	  std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	  if(faceNodes[0]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[0]->getLocalID() != faceNodes1[1]->getLocalID()){
	    found=true;
	    thirdNode =faceNodes[0];
	  }
	if(faceNodes[1]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[1]->getLocalID() != faceNodes1[1]->getLocalID()){
	  found=true;
	  thirdNode =faceNodes[1];
	}
	}
      }
      for (CFuint iDim=0; iDim<nbDims ;++iDim){
	v1[iDim] =(*faceNodes1[0])[iDim]-(*faceNodes1[1])[iDim];
	v2[iDim] =(*faceNodes1[0])[iDim]-(*thirdNode)[iDim];
	v3[iDim]=(*faceNodes1[1])[iDim]-(*thirdNode)[iDim];
      }
      MathTools::MathFunctions::crossProd(v1, v2,res);
      
      CFreal area =res.norm2();
      CFreal perim = v1.norm2()+v2.norm2()+v3.norm2();
      CFreal radius =area/perim;
      m_mapNodeRadius.insert(faceNodes1[0]->getLocalID(),radius);
      m_mapNodeNState.insert(faceNodes1[0]->getLocalID() ,v1state);
      
      m_mapNodeRadius.insert(faceNodes1[1]->getLocalID(),radius);
      m_mapNodeNState.insert(faceNodes1[1]->getLocalID(),v1state);
      
      m_mapNodeRadius.insert(thirdNode->getLocalID(),radius);
      m_mapNodeNState.insert(thirdNode->getLocalID(),v1state);
      
      m_geoBuilder.releaseGE();
      m_mapNodeRadius.sortKeys();
      m_mapNodeNState.sortKeys(); 
    }
    m_mapNodeRadius1=m_mapNodeRadius;
    m_mapNodeNState1=m_mapNodeNState;
  }
  
  
  
  /////////////////////////
  //2D quadrilateral mesh//
  /////////////////////////
  
  // Connectivity information 2D quadrilateral
  if (nbElemTypes==4 && nbDims==2){
    CFuint nbPairsNodeNode = 200000; 
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeNode(nbPairsNodeNode);
  
  std::multimap<CFuint, CFuint>  mapNodeNode;
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  CFuint coupleDone [nbPairsNodeNode][2];
  for (CFuint k= 0 ; k<nbPairsNodeNode ; ++k){
    coupleDone[k][0]=0;
    coupleDone[k][1]=0;
  }
  CFuint i = 0;
  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size(); 

      for (CFuint iFace=0; iFace<nbFaces; ++iFace){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	  CFuint nodeIDinF1 = faceNodes[0]->getLocalID();
	  CFuint nodeIDinF2 = faceNodes[1]->getLocalID();

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

      m_geoBuilder.releaseGE();
  }

  m_mapNodeNode1=m_mapNodeNode;
  
  // End of 2D quadrilateral connectivity information

  ////////////////////////////////////
  // 2D quads MQI based on skewness //
  ////////////////////////////////////
  if (nbElemTypes==4 && nbDims==2 && m_MQIvalue==4){
  CFuint nbPairsNodeSkew =300000 ; 
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  typedef CFMultiMap<CFuint, CFreal> MapNodeNState;
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNState(nbPairsNodeSkew);
  typedef CFMultiMap<CFuint, CFreal> MapNodeSkew;
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeSkew(nbPairsNodeSkew);
  typedef MapNodeNode::MapIterator mapIt;
  typedef MapNodeNode::MapIterator mapItN;

  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    Framework::Node* firstNode;
    Framework::Node* secondNode;
    Framework::Node* thirdNode;
    Framework::Node* fourthNode;

    const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size();
    //
    //
    // std::cout<<facesInCell<<endl;
    //
    std::vector<Framework::Node*>& faceNodes12 = *facesInCell[0]->getNodes();
    firstNode = faceNodes12[0];
    secondNode = faceNodes12[1];
    bool foundFN =false;
    bool foundSN= false;
    std::pair<mapIt,mapIt > itFirstNode=m_mapNodeNode1.find(firstNode->getLocalID(), foundFN);
    std::pair<mapItN,mapItN > itSecondNode=m_mapNodeNode1.find(secondNode->getLocalID(), foundSN);
    cf_assert(foundFN);
    cf_assert(foundSN);
    for (CFuint iFace=1; iFace<nbFaces; ++iFace){
      std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
      const CFuint nbNodesinF = faceNodes.size();
      for (mapIt itFN = itFirstNode.first; itFN != itFirstNode.second; ++itFN) {
	for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	  if(itFN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	    fourthNode =  faceNodes[iNode];
	  }
	}
      }      
      for (mapItN itSN = itSecondNode.first; itSN != itSecondNode.second; ++itSN) {
	for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	  if(itSN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	    thirdNode =  faceNodes[iNode];
	  }
	}
      }
    }

  //First order simulation v1state = user-defined
  //  variable to be putted w.r.t the monitor variable
  //  initial state  
    CFreal v1state=1000.;

  ////////////////////////////////////////////////////
  //Starting from second order v1state = NodalStates//
  ////////////////////////////////////////////////////
    /*   if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
      v1state= (*states[iCell])[m_monitorVarID];
    }
    else{
      cf_assert(m_monitorPhysVarID < m_pdata.size());
      m_state->copyData(*states[iCell]);
      getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
      v1state  = m_pdata[m_monitorPhysVarID];
      }*/
    //
    //End of second order
    //
    
     
    CFreal skew =computeSkewness2dQuads(firstNode,secondNode,thirdNode,fourthNode);
    m_mapNodeSkew.insert(firstNode->getLocalID(),skew);
    m_mapNodeSkew.insert(secondNode->getLocalID(),skew);
    m_mapNodeSkew.insert(thirdNode->getLocalID(),skew);
    m_mapNodeSkew.insert(fourthNode->getLocalID(),skew);

    m_mapNodeNState.insert(firstNode->getLocalID(),v1state);
    m_mapNodeNState.insert(secondNode->getLocalID(),v1state);
    m_mapNodeNState.insert(thirdNode->getLocalID(),v1state);
    m_mapNodeNState.insert(fourthNode->getLocalID(),v1state);

    m_mapNodeSkew.sortKeys();
    m_mapNodeNState.sortKeys();
    m_geoBuilder.releaseGE();

  }


    m_mapNodeSkew1=m_mapNodeSkew;
    m_mapNodeNState1=m_mapNodeNState;
  }
 
  ////////////////////////////////////////
  // 2D quads MQI based on Aspect ratio //
  ////////////////////////////////////////

  if (nbElemTypes==4 && nbDims==2 &&  m_MQIvalue==3){
  CFuint nbPairsNodeAR =300000 ; 
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  typedef CFMultiMap<CFuint, CFreal> MapNodeNState;
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNState(nbPairsNodeAR);
  typedef CFMultiMap<CFuint, CFreal> MapNode;

  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeAR(nbPairsNodeAR);
  typedef MapNodeNode::MapIterator mapIt;
  typedef MapNodeNode::MapIterator mapItN;

  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    Framework::Node* firstNode;
    Framework::Node* secondNode;
    Framework::Node* thirdNode;
    Framework::Node* fourthNode;

    const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size(); 
    std::vector<Framework::Node*>& faceNodes12 = *facesInCell[0]->getNodes();
    firstNode = faceNodes12[0];
    secondNode = faceNodes12[1];
    bool foundFN =false;
    bool foundSN= false;
    std::pair<mapIt,mapIt > itFirstNode=m_mapNodeNode1.find(firstNode->getLocalID(), foundFN);
    std::pair<mapItN,mapItN > itSecondNode=m_mapNodeNode1.find(secondNode->getLocalID(), foundSN);
    cf_assert(foundFN);
    cf_assert(foundSN);
    for (CFuint iFace=1; iFace<nbFaces; ++iFace){
      std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
      const CFuint nbNodesinF = faceNodes.size();
      for (mapIt itFN = itFirstNode.first; itFN != itFirstNode.second; ++itFN) {
	for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	  if(itFN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	    fourthNode =  faceNodes[iNode];
	  }
	}
      }      
      for (mapItN itSN = itSecondNode.first; itSN != itSecondNode.second; ++itSN) {
	for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	  if(itSN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	    thirdNode =  faceNodes[iNode];
	  }
	}
      }
    }
     
  //First order simulation v1state = user-defined
  //  variable to be putted w.r.t the monitor variable
  //  initial state  , It is an example of a particular test case
    CFreal v1state =0.00515086;

  ////////////////////////////////////////////////////
  //Starting from second order v1state = NodalStates//
  ////////////////////////////////////////////////////
    /*   if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
      v1state= (*states[iCell])[m_monitorVarID];
    }
    else{
      cf_assert(m_monitorPhysVarID < m_pdata.size());
      m_state->copyData(*states[iCell]);
      getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
      v1state  = m_pdata[m_monitorPhysVarID];
      }*/
    //
    //End of second order
    //

    CFreal AR =computeAspectRatio2dQuads(firstNode,secondNode,thirdNode,fourthNode);

    m_mapNodeAR.insert(firstNode->getLocalID(),AR);
    m_mapNodeAR.insert(secondNode->getLocalID(),AR);
    m_mapNodeAR.insert(thirdNode->getLocalID(),AR);
    m_mapNodeAR.insert(fourthNode->getLocalID(),AR);

    m_mapNodeNState.insert(firstNode->getLocalID(),v1state);
    m_mapNodeNState.insert(secondNode->getLocalID(),v1state);
    m_mapNodeNState.insert(thirdNode->getLocalID(),v1state);
    m_mapNodeNState.insert(fourthNode->getLocalID(),v1state);

    m_mapNodeAR.sortKeys();
    m_mapNodeNState.sortKeys();
    m_geoBuilder.releaseGE();

    }
    m_mapNodeAR1=m_mapNodeAR;
    m_mapNodeNState1=m_mapNodeNState;
  }
  
  }

    ////////////////////////////////////////////
    //Connectivity information 3D teterahedral//
    ////////////////////////////////////////////
  if (nbElemTypes==4 && nbDims==3 ){
   Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
 
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  CFuint nbPairsNodeNode1 =300000; 
  typedef CFMultiMap<CFuint, CFuint> MapCellNode;
  Common::CFMultiMap<CFuint,CFuint>  m_mapCellNode(nbPairsNodeNode1);

  typedef CFMultiMap<CFuint, CFreal > MapNodeNState;
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNState(nbPairsNodeNode1);

  typedef CFMultiMap<CFuint, CFreal > MapNodeRS;
  Common::CFMultiMap<CFuint,CFreal>  m_mapNodeTS(nbPairsNodeNode1);
  typedef MapCellNode::MapIterator mapItc;

  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
    for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
      CFuint nodeIDinC=cells->getNodeID(iCell, iNodeC);
      m_mapCellNode.insert(iCell,nodeIDinC);
      }
    m_mapCellNode.sortKeys();

    m_geoBuilder.releaseGE();
  }
  m_mapCellNode1=m_mapCellNode;
  //end of 3D connectivity information

  //////////////////////////////////////////
  //MQI based on radius of inserted sphere//
  //////////////////////////////////////////

  if (nbElemTypes==4 && nbDims==3 && m_MQIvalue==5){
   for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    std::vector<CFuint> nodeIDs;
    nodeIDs.clear();
    CFreal facesArea = 0.;
    CFreal volume =0.;
    CFreal radius = 0.;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    bool foundC = false;
    std::pair<mapItc,mapItc > iteC=m_mapCellNode1.find(iCell, foundC);
    cf_assert(foundC);
    for (mapItc it2 = iteC.first; it2 != iteC.second; ++it2){
      nodeIDs.push_back(it2->second);
    }
    volume=ComputeTvolume(nodeIDs[0],nodeIDs[1],nodeIDs[2],nodeIDs[3]);
    const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size(); 
    for (CFuint iFace=1; iFace<nbFaces; ++iFace){
      std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
      Framework::Node* node1 = faceNodes[0];
      Framework::Node* node2 = faceNodes[1];
      Framework::Node* node3 = faceNodes[2];
      facesArea+= ComputeTFacesurface(node1,node2,node3);
    }
    radius = 3.*volume/facesArea;

     
  //First order simulation v1state = user-defined
  //  variable to be putted w.r.t the monitor variable
  //  initial state  
    CFreal v1state =1000.;
  ////////////////////////////////////////////////////
  //Starting from second order v1state = NodalStates//
  ////////////////////////////////////////////////////
    /*   if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
      v1state= (*states[iCell])[m_monitorVarID];
    }
    else{
      cf_assert(m_monitorPhysVarID < m_pdata.size());
      m_state->copyData(*states[iCell]);
      getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
      v1state  = m_pdata[m_monitorPhysVarID];
      }*/
    //
    //End of second order
    //

     for(CFuint i=0; i<nodeIDs.size() ; ++i){
      m_mapNodeTS.insert(nodeIDs[i],radius);
      m_mapNodeNState.insert(nodeIDs[i],v1state);
    }
    m_mapNodeTS.sortKeys();
    m_mapNodeNState.sortKeys();

    m_geoBuilder.releaseGE();

 }
    m_mapNodeTS1=m_mapNodeTS;
    m_mapNodeNState1=m_mapNodeNState;

  }
}

   /////////////////////////
  //  3D hexahedral mesh //
  /////////////////////////
  
  // Connectivity information 3D hexahedral
  if (nbElemTypes==8 && nbDims==3){
    time_t tstart, tend; 
    tstart = time(0);
    Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
    const CFuint nbNodes = nodes.size();
    CFuint nbPairsNodeNode =  6*nbNodes;
    typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
    Common::CFMultiMap<CFuint,CFuint>  m_mapNodeNode(5000000);
    typedef MapNodeNode::MapIterator mapNodeIt;
    
    Common::SafePtr<Framework::TopologicalRegionSet> cells = 
      Framework::MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    typedef CFMultiMap<CFuint, CFuint> MapCellNode;
    Common::CFMultiMap<CFuint,CFuint>  m_mapCellNode(5000000);
    typedef MapCellNode::MapIterator mapCellIt;
    
    typedef CFMultiMap<CFuint, CFuint> MapNodeCell;
    Common::CFMultiMap<CFuint,CFuint>  m_mapNodeCell(5000000);
    typedef MapNodeCell::MapIterator mapNodeIt;
    
    typedef CFMultiMap<CFuint, CFuint > MapFaceNode;
    Common::CFMultiMap<CFuint,CFuint>  m_mapFaceNode(5000000);
    typedef MapFaceNode::MapIterator mapFaceIt;
    
    std::vector<CFuint> faceIDs;
    std::vector<CFuint> nodesInFaceIDs;
    std::vector<CFuint> cellsInNodeIDs;
    std::vector<CFuint> cellsInNeighborIDs;	
    
    CFuint coupleDone [nbPairsNodeNode][2];
    for (CFuint k= 0 ; k<nbPairsNodeNode ; ++k){
      coupleDone[k][0]=0;
      coupleDone[k][1]=0;
    }
    
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
     geoData.trs = cells;
     
     
     for (CFuint iCell=0; iCell<nbCells; ++iCell){
       geoData.idx = iCell;
       Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
       const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
       for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
	 CFuint nodeIDinC=cells->getNodeID(iCell, iNodeC);
	 m_mapCellNode.insert(iCell,nodeIDinC);
	 m_mapNodeCell.insert(nodeIDinC, iCell);
       }
       m_mapNodeCell.sortKeys();
       
       m_geoBuilder.releaseGE();
     }
     CFuint i = 0;
     for (CFuint iCell=0; iCell<nbCells; ++iCell){
       geoData.idx = iCell;
       Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
       const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
       const CFuint nbFaces = facesInCell.size(); 
       for (CFuint iFace=0; iFace<nbFaces; ++iFace){
	 CFuint faceID = facesInCell[iFace]->getID();
	 
	 bool faceFound = false;
	 std::pair<mapFaceIt,mapFaceIt > ite=m_mapFaceNode.find(faceID, faceFound);
	 if (!faceFound){
	   bool isGhost = facesInCell[iFace]->getState(1)->isGhost();
	   std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	   for (CFuint iNodeC = 0; iNodeC < faceNodes.size(); ++iNodeC){
	     m_mapFaceNode.insert(faceID,faceNodes[iNodeC]->getLocalID());
	   }  
	   for (CFuint k =0; k<faceNodes.size()-1; ++k){
	     cellsInNodeIDs.clear();
	     CFuint nodeID = faceNodes[k]->getLocalID();
	     for (mapNodeIt it = m_mapNodeCell.begin(); it!= m_mapNodeCell.end(); ++it){
	       if (it->first == nodeID){
		 cellsInNodeIDs.push_back(it->second);
	       }
	     }
	     for(CFuint l= k+1; l<faceNodes.size(); ++l){
	       cellsInNeighborIDs.clear();
	       CFuint neighborNodeID = faceNodes[l]->getLocalID();
	       CFuint cellsInCommon = 0;
	       for (mapNodeIt it = m_mapNodeCell.begin(); it!= m_mapNodeCell.end(); ++it){
		 if (it->first == neighborNodeID){
		   cellsInNeighborIDs.push_back(it->second);
		 }
	       }
	       for (CFuint m = 0; m<cellsInNodeIDs.size(); ++m){
		 for (CFuint n = 0; n<cellsInNeighborIDs.size(); ++n){	
		   if (cellsInNodeIDs[m] == cellsInNeighborIDs[n]){
		     cellsInCommon +=1;
		   }
		 }
	       }
	       if ((cellsInCommon == 4) || (cellsInCommon == 2 && isGhost) || (isGhost && cellsInCommon == 1 && ((cellsInNeighborIDs.size() == 1 && cellsInNodeIDs.size() == 2) ||  (cellsInNeighborIDs.size() == 2 &&   cellsInNodeIDs.size() == 1))) || (isGhost && cellsInCommon == 1 && cellsInNeighborIDs.size() == 2 && cellsInNodeIDs.size() == 2)){
		 bool done1 = false;
		 bool done2 = false;
		 CFuint nodeIDinF1 = nodeID;
		 CFuint nodeIDinF2 = neighborNodeID;
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
	     }
	   }
	   
	       }   
       }
       m_mapCellNode.sortKeys();
       m_mapNodeNode.sortKeys();
       m_mapFaceNode.sortKeys();
       m_geoBuilder.releaseGE();
     }
     m_mapNodeNode1 = m_mapNodeNode;
     m_mapCellNode1 = m_mapCellNode;
     m_mapNodeCell1 = m_mapNodeCell;
     tend = time(0); 
     cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;	
  }
  
// End of 3D hexahedral connectivity information
  if(m_interpolateState){
    CFLog(INFO, "Interpolation Activated ---- Save Old Mesh properties\n");
    saveOldMeshProperties();
  }	
}
      
      //////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::unsetup()
{
  CFAUTOTRACE;

  m_geoBuilder.unsetup();
  m_faceTRSBuilder.unsetup();
  
  Framework::DataProcessingCom::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////

 void MeshFittingAlgorithm::createNodalConnectivity()
{ 
  m_edgeGraph.setNodeDataSocket(socket_nodes);
  m_edgeGraph.computeConnectivity(); 
  m_edgeGraphN.setNodeDataSocket(socket_nodes);
  m_edgeGraphN.computeConnectivity(); 
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
      if (facesInCell[iFace]->getState(1)->isGhost() ) {
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
        CFLog(DEBUG_MED,"averageNormal  "<<averageNormal<< "/n");
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
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
 
  typedef CFMultiMap<CFuint, CFuint> MapNodeCell;
  Common::CFMultiMap<CFuint,CFuint>  m_mapNodeCell(5000000); // AL: what is this??? please compute some estimation of size which is more automatic 
  typedef MapNodeCell::MapIterator mapIt;
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
 	nodeDistance[nodes[iNode]->getLocalID()] = 0.;
  }
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
	CFuint adjacentCells = 0;
	for (mapIt it = m_mapNodeCell1.begin(); it!= m_mapNodeCell1.end(); ++it){
		if (it->first==nodes[iNode]->getLocalID()){
			nodeDistance[nodes[iNode]->getLocalID()]+= wallDistance[it->second];
			adjacentCells +=1;
		}
         }
	nodeDistance[nodes[iNode]->getLocalID()] = nodeDistance[nodes[iNode]->getLocalID()]/adjacentCells;
   }
	
			
  
  movingIter+=1;

  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::execute() => start \n");
  resizeSystemSolverToNodalData();
  computeSpringTruncationData(); 
  solveLinearSystem();
    Common::SafePtr<Framework::TopologicalRegionSet> cells =
      Framework::MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbElemTypes = cells->getNbNodesInGeo(0);

  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  Framework::DataHandle<CFreal>relativeError  = socket_relativeError.getDataHandle();

  if(  m_tolerance!=0. &&  nbElemTypes == 4 && nbDims == 2 ){
    Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();  
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
    typedef MapNodeNode::MapIterator mapIt;
    typedef MapNodeNode::MapIterator mapItN;
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    typedef CFMultiMap<CFuint, CFreal> MapNodeARF;
    // FB:hard coded, need to be revised 
    CFuint nbPairsNodeAR =200000 ; 
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeARF(nbPairsNodeAR);
    geoData.trs = cells;
    CFuint nbPairsNodeNode = 300000;
    typedef CFMultiMap<CFuint, CFreal> MapNodeRadiusF;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeRadiusF(nbPairsNodeNode);
    typedef CFMultiMap<CFuint, CFreal> MapNodeARF1;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeARF1(nbPairsNodeAR);
    ////////////////////////////////
    //RSI for quadrilateral meshes//
    //RSI computations based on AR//
    //RSI for time n              //
    ////////////////////////////////
    
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
      geoData.idx = iCell;
      Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
      Framework::Node * firstNode;
      Framework::Node * secondNode;
      Framework::Node * thirdNode;
      Framework::Node * fourthNode;
      const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
      const CFuint nbFaces = facesInCell.size(); 
      std::vector<Framework::Node*>& faceNodes12 = *facesInCell[0]->getNodes();
      firstNode = faceNodes12[0];
      secondNode = faceNodes12[1];
      bool foundFN =false;
      bool foundSN =false;
      std::pair<mapIt,mapIt > itFirstNode=m_mapNodeNode1.find(firstNode->getLocalID(), foundFN);
      std::pair<mapItN,mapItN > itSecondNode=m_mapNodeNode1.find(secondNode->getLocalID(), foundSN);
      cf_assert(foundFN);
      cf_assert(foundSN);
      for (CFuint iFace=1; iFace<nbFaces; ++iFace){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	const CFuint nbNodesinF = faceNodes.size();
	for (mapIt itFN = itFirstNode.first; itFN != itFirstNode.second; ++itFN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itFN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	      fourthNode =  faceNodes[iNode];
	    }
	  }
	}      
	for (mapItN itSN = itSecondNode.first; itSN != itSecondNode.second; ++itSN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itSN->second == faceNodes[iNode]->getLocalID()  && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID() ){
	      thirdNode =  faceNodes[iNode];
	    }
	  }
	}
      }
      CFreal AR =computeAspectRatio2dQuads(firstNode,secondNode,thirdNode,fourthNode);
      m_mapNodeARF.insert(firstNode->getLocalID(),AR);
      m_mapNodeARF.insert(secondNode->getLocalID(),AR);
      m_mapNodeARF.insert(thirdNode->getLocalID(),AR);
      m_mapNodeARF.insert(fourthNode->getLocalID(),AR);
      
      m_mapNodeARF.sortKeys();

      m_geoBuilder.releaseGE();
      
    }

    // mesh node repositionning
    updateNodePositions();
    //

    ////////////////////////////////
    //RSI for quadrilateral mesh  //
    //RSI computations based on AR//
    //RSI for time n+1            //
    ////////////////////////////////
    geoData.trs = cells;
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
      geoData.idx = iCell;
      Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
      Framework::Node * firstNode;
      Framework::Node * secondNode;
      Framework::Node * thirdNode;
      Framework::Node * fourthNode;
      
      const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
      const CFuint nbFaces = facesInCell.size(); 
      std::vector<Framework::Node*>& faceNodes12 = *facesInCell[0]->getNodes();
      firstNode = faceNodes12[0];
      secondNode = faceNodes12[1];
      bool foundFN =false;
      bool foundSN =false;
      std::pair<mapIt,mapIt > itFirstNode=m_mapNodeNode1.find(firstNode->getLocalID(), foundFN);
      std::pair<mapItN,mapItN > itSecondNode=m_mapNodeNode1.find(secondNode->getLocalID(), foundSN);
      cf_assert(foundFN);
      cf_assert(foundSN);
      for (CFuint iFace=1; iFace<nbFaces; ++iFace){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	const CFuint nbNodesinF = faceNodes.size();
	for (mapIt itFN = itFirstNode.first; itFN != itFirstNode.second; ++itFN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itFN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	    fourthNode =  faceNodes[iNode];
	    }
	  }
	}      
	for (mapItN itSN = itSecondNode.first; itSN != itSecondNode.second; ++itSN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itSN->second == faceNodes[iNode]->getLocalID()  && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID() ){
	      thirdNode =  faceNodes[iNode];
	    }
	  }
	}
      }
      
      
      CFreal AR =computeAspectRatio2dQuads(firstNode,secondNode,thirdNode,fourthNode);
      m_mapNodeARF1.insert(firstNode->getLocalID(),AR);
      m_mapNodeARF1.insert(secondNode->getLocalID(),AR);
      m_mapNodeARF1.insert(thirdNode->getLocalID(),AR);
      m_mapNodeARF1.insert(fourthNode->getLocalID(),AR);
      
      m_mapNodeARF1.sortKeys();
      
      m_geoBuilder.releaseGE();
      
    }
    typedef MapNodeARF1::MapIterator mapItSkf;
    typedef MapNodeARF::MapIterator mapItSki;
    CFreal Rn=0;
    CFreal Rn1=0;
    CFreal dispup=0.;
    CFreal dispdown=0.;
    CFreal dispmoyup=0.;
    CFreal dispmoydown=0.;
    CFreal a=0.;
    CFreal b=0.;
    for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
      if (nodes[iNode]->isParUpdatable() ){  //&& isNodeLocked(nodes[iNode]) == false){
	bool foundRi = false;
	bool foundRf = false;
	std::pair<mapItSkf,mapItSkf > itrf=m_mapNodeARF1.find(nodes[iNode]->getLocalID(), foundRf);
	std::pair<mapItSki,mapItSki > itri=m_mapNodeARF.find(nodes[iNode]->getLocalID(), foundRi);
	cf_assert(foundRf);
	cf_assert(foundRi);
	for (mapItSki it = itri.first; it != itri.second; ++it){
	  Rn +=it->second;
	}
	for (mapItSkf it = itrf.first; it != itrf.second; ++it){
	  Rn1 +=it->second;
	}
	CFreal disp=std::abs(Rn1-Rn)/Rn;
	if(disp>0.01){
	  dispup+=disp;
	  a+=1.;
	}
	if(disp <= 0.01){
	  dispdown+=disp;
	  b+=1.;
	}
      }
    }
    
    if(a > 0.){
      dispmoyup=dispup/a;
    }
    if(b> 0.){
      dispmoydown=dispdown/b;
    }
    
    CFreal x= std::sqrt(a)*dispmoyup;
    CFreal erreur= (pow(a,x)*std::sqrt(dispmoyup)+b*dispmoydown)/nodes.size();
    // CFreal epsilonFinal=0.;
    const std::string nsp = this->getMethodData().getNamespace();
    MPI_Comm communicator = Common::PE::GetPE().GetCommunicator(nsp); 
    
    MPI_Allreduce(&erreur,&m_erreurG, 1,MPI_DOUBLE,MPI_MAX, communicator);
    
    for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
      relativeError[nodes[iNode]->getLocalID()] =m_erreurG;	
    }
  }
  

  ////////////////////////////////////
  //RSI for triangular  meshes      //
  //RSI computations based on radius//
  //RSI for time n                  //
  ////////////////////////////////////

  if( m_tolerance!=0. &&  nbElemTypes == 3 && nbDims == 2 ){
    Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();  
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
    typedef MapNodeNode::MapIterator mapIt;
    typedef MapNodeNode::MapIterator mapItN;
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    geoData.trs = cells;
    CFuint nbPairsNodeNode = 300000;
    typedef CFMultiMap<CFuint, CFreal> MapNodeRadiusF;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeRadiusF(nbPairsNodeNode);
    typedef CFMultiMap<CFuint, CFreal> MapNodeRadiusF1;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeRadiusF1(nbPairsNodeNode);
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    Framework::Node * thirdNode;
    const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size(); 
    RealVector res(3);res=0.;
    std::vector<Framework::Node*>& faceNodes1 = *facesInCell[0]->getNodes();
    bool found=false;
    RealVector v1(3); v1=0.;
    RealVector v2(3); v2=0.;
    RealVector v3(3); v3=0.;
    //CFreal v1state;
    for (CFuint iFace=1; iFace<nbFaces; ++iFace){
      while(found==false){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	if(faceNodes[0]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[0]->getLocalID() != faceNodes1[1]->getLocalID()){
	  found=true;
	  thirdNode =faceNodes[0];
	}
	if(faceNodes[1]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[1]->getLocalID() != faceNodes1[1]->getLocalID()){
	  found=true;
	  thirdNode =faceNodes[1];
	}
      }
    }
    for (CFuint iDim=0; iDim<nbDims ;++iDim){
      v1[iDim] =(*faceNodes1[0])[iDim]-(*faceNodes1[1])[iDim];
      v2[iDim] =(*faceNodes1[0])[iDim]-(*thirdNode)[iDim];
      v3[iDim]=(*faceNodes1[1])[iDim]-(*thirdNode)[iDim];
    }
    MathTools::MathFunctions::crossProd(v1, v2,res);
    CFreal area =res.norm2();
    CFreal perim = v1.norm2()+v2.norm2()+v3.norm2();
    CFreal radius =area/perim;

    m_mapNodeRadiusF.insert(faceNodes1[0]->getLocalID(),radius);
    m_mapNodeRadiusF.insert(faceNodes1[1]->getLocalID(),radius);
    m_mapNodeRadiusF.insert(thirdNode->getLocalID(),radius);

    m_geoBuilder.releaseGE();
    m_mapNodeRadiusF.sortKeys();
  
    }

    // mesh node repositionning
    updateNodePositions();
    //
  ////////////////////////////////////
  //RSI for triangular  meshes      //
  //RSI computations based on radius//
  //RSI for time n+1                //
  ////////////////////////////////////

  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
    Framework::Node * thirdNode;
    const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
    const CFuint nbFaces = facesInCell.size(); 
    RealVector res(3);res=0.;
    std::vector<Framework::Node*>& faceNodes1 = *facesInCell[0]->getNodes();
    bool found=false;
    RealVector v1(3); v1=0.;
    RealVector v2(3); v2=0.;
    RealVector v3(3); v3=0.;
    // CFreal v1state;
    for (CFuint iFace=1; iFace<nbFaces; ++iFace){
      while(found==false){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	if(faceNodes[0]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[0]->getLocalID() != faceNodes1[1]->getLocalID()){
	  found=true;
	  thirdNode =faceNodes[0];
	}
	if(faceNodes[1]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[1]->getLocalID() != faceNodes1[1]->getLocalID()){
	  found=true;
	  thirdNode =faceNodes[1];
	}
      }
    }
    for (CFuint iDim=0; iDim<nbDims ;++iDim){
      v1[iDim] =(*faceNodes1[0])[iDim]-(*faceNodes1[1])[iDim];
      v2[iDim] =(*faceNodes1[0])[iDim]-(*thirdNode)[iDim];
      v3[iDim]=(*faceNodes1[1])[iDim]-(*thirdNode)[iDim];
    }
    MathTools::MathFunctions::crossProd(v1, v2,res);
    
    CFreal area =res.norm2();
    CFreal perim = v1.norm2()+v2.norm2()+v3.norm2();
    CFreal radius =area/perim;


    m_mapNodeRadiusF1.insert(faceNodes1[0]->getLocalID(),radius);
    m_mapNodeRadiusF1.insert(faceNodes1[1]->getLocalID(),radius);
    m_mapNodeRadiusF1.insert(thirdNode->getLocalID(),radius);

    m_geoBuilder.releaseGE();
    m_mapNodeRadiusF1.sortKeys();
  
    }

  typedef MapNodeRadiusF1::MapIterator mapItSkf;
  typedef MapNodeRadiusF::MapIterator mapItSki;
  CFreal Rn=0.;
  CFreal Rn1=0;
  CFreal dispup=0.;
  CFreal dispdown=0.;
  CFreal dispmoyup=0.;
  CFreal dispmoydown=0.;
  CFreal a=0.;
  CFreal b=0.;
   for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
     if (nodes[iNode]->isParUpdatable() ){
       bool foundRi = false;
       bool foundRf = false;
       std::pair<mapItSkf,mapItSkf > itrf=m_mapNodeRadiusF1.find(nodes[iNode]->getLocalID(), foundRf);
       std::pair<mapItSki,mapItSki > itri=m_mapNodeRadiusF.find(nodes[iNode]->getLocalID(), foundRi);
       cf_assert(foundRf);
       cf_assert(foundRi);
       for (mapItSki it = itri.first; it != itri.second; ++it){
	 Rn +=it->second;
       }
       for (mapItSkf it = itrf.first; it != itrf.second; ++it){
	 Rn1 +=it->second;
       }
       CFreal disp=std::abs(Rn1-Rn)/Rn;
       CFreal ref= m_tolerance/100.;
       if(disp>ref){
	 dispup+=disp;
	 a+=1.;
       }
       if(disp <= ref){
	 dispdown+=disp;
	 b+=1.;
     }
   }
   }
   if(a > 0.){
   dispmoyup=dispup/a;
   }
   if(b> 0.){
   dispmoydown=dispdown/b;
   }
   CFreal x= std::sqrt(a)*dispmoyup;
   CFreal erreur= (pow(a,x)*std::sqrt(dispmoyup)+b*dispmoydown)/nodes.size();
   CFreal epsilonFinal=0.;
   const std::string nsp = this->getMethodData().getNamespace();
   // const int nbProcesses = Common::PE::GetPE().GetProcessorCount(nsp);
   const int processRank = Common::PE::GetPE().GetRank(nsp);
   MPI_Comm communicator = Common::PE::GetPE().GetCommunicator(nsp); 
   
   MPI_Allreduce(&erreur,&m_erreurG, 1,MPI_DOUBLE,MPI_MAX, communicator);

   /* if(processRank==0){
     epsilonFinal=m_erreurG;
     std::ofstream myFile ; 
     myFile.open("Quad.txt");
     myFile << m_erreurG << "      " << SubSystemStatusStack::getActive()->getNbIter() << "\n";
     myFile.close();
     }*/
   for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
     relativeError[nodes[iNode]->getLocalID()]=m_erreurG;
   }
  }
  const std::string nsp = this->getMethodData().getNamespace();
  const int nbProcesses = Common::PE::GetPE().GetProcessorCount(nsp);
  const int processRank = Common::PE::GetPE().GetRank(nsp);
  MPI_Comm communicator = Common::PE::GetPE().GetCommunicator(nsp);
  //
  if(processRank==0){
  
  std::ofstream myFile ;
  myFile.open("Quad.txt", ios::app);
  myFile << m_erreurG << "  " << SubSystemStatusStack::getActive()->getNbIter() << "\n";
  myFile.close();
  
  }
  
  ///////////////////
  //RSI deactivated//
  ///////////////////

  if(m_tolerance == 0.){
    updateNodePositions(); 
  }
  

  resizeSystemSolverToStateData();
  triggerRecomputeMeshData();
  // stateInterpolator();
  /*if(m_interpolateState){ 
      interpolateNewMeshProperties();     
      CFLog(INFO, "Interpolation Step ==> end\n");
      }*/
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
  CFreal MeshFittingAlgorithm::computeAspectRatio2dQuads(const  Framework::Node* const  a,const  Framework::Node* const  b,
							 const  Framework::Node* const  c, const  Framework::Node* const  d){
            
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector vectorDC(3); vectorDC = 0.;
  RealVector vectorDA(3); vectorDA = 0.;
  RealVector vectorBC(3); vectorBC = 0.;
  RealVector vectorBA(3); vectorBA = 0.;
  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    vectorDC[iDim]= (*d)[XX+iDim]-(*c)[XX+iDim];
    vectorDA[iDim]= (*d)[XX+iDim]-(*a)[XX+iDim];
    vectorBC[iDim]= (*b)[XX+iDim]-(*c)[XX+iDim];
    vectorBA[iDim]= (*b)[XX+iDim]-(*a)[XX+iDim];
  }
  
  CFreal min1=vectorBA.norm2();
  CFreal AR =vectorBC.norm2()/vectorBA.norm2();
  if(min1>vectorBC.norm2()){
    AR =vectorBA.norm2()/vectorBC.norm2();
  }
  return AR;
}
//////////////////////////////////////////////////////////////////////////////
     
   CFreal MeshFittingAlgorithm::computeSkewness2dQuads(const  Framework::Node* const  a1,const  Framework::Node* const  b1,
						       const  Framework::Node* const  c1, const  Framework::Node* const  d1){

   Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();   
   const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
   RealVector vectorDC(nbDims); vectorDC = 0.;
   RealVector vectorDA(nbDims); vectorDA = 0.;
   RealVector vectorAD(nbDims); vectorAD = 0.;
   RealVector vectorAB(nbDims); vectorAB = 0.;
   RealVector vectorCD(nbDims); vectorCD = 0.;
   RealVector vectorBA(nbDims); vectorBA = 0.;
   RealVector vectorBC(nbDims); vectorBC = 0.;
   RealVector vectorCB(nbDims); vectorCB = 0.;

   CFreal angle1;
   CFreal angle2;
   CFreal angleC;
   CFreal angleB;


   CFreal skew;
   CFreal skew2;
   for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
     vectorDC[iDim]= (*c1)[XX+iDim]-(*d1)[XX+iDim];
     vectorDA[iDim]= (*a1)[XX+iDim]-(*d1)[XX+iDim];

     vectorAB[iDim]= (*b1)[XX+iDim]-(*a1)[XX+iDim];
     vectorAD[iDim]= (*d1)[XX+iDim]-(*a1)[XX+iDim];

     vectorCD[iDim]= -vectorDC[XX+iDim];
     vectorBA[iDim]= -vectorAB[XX+iDim];
   
     vectorBC[iDim]= (*c1)[XX+iDim]-(*b1)[XX+iDim];
     vectorCB[iDim]= -vectorBC[XX+iDim];

   }
   CFreal cosADC = MathTools::MathFunctions::innerProd(vectorDC,vectorDA)/(vectorDC.norm2()*vectorDA.norm2());
     angle1 = acos(cosADC);
    
   CFreal cosDCB = MathTools::MathFunctions::innerProd(vectorCD,vectorCB)/(vectorCD.norm2()*vectorCB.norm2());
     angleC=acos(cosDCB);

   CFreal cosABC = MathTools::MathFunctions::innerProd(vectorBA,vectorBC)/(vectorBA.norm2()*vectorBC.norm2());
     angleB = acos(cosABC);
   CFreal cosDAB = MathTools::MathFunctions::innerProd(vectorAB,vectorAD)/(vectorAB.norm2()*vectorAD.norm2());
   angle2 = acos(cosDAB);
   
   CFreal pi2 = 3.1415/2.0;
    if(angle2>angle1){
     skew =std:: max((angle2-pi2)/pi2  ,  (pi2-angle1)/pi2 );
   }
   else{
     skew =std:: max((angle1-pi2)/pi2  ,  (pi2-angle2)/pi2 );
     }

    if(angleC>angleB){
     skew2 =std:: max((angleC-pi2)/pi2  ,  (pi2-angleB)/pi2 );
   }
   else{
     skew2 =std:: max((angleB-pi2)/pi2  ,  (pi2-angleC)/pi2 );
     }

   CFreal finalSkew = std::max( std::abs(skew), std::abs(skew2)); 

   return finalSkew;
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

  CFreal MeshFittingAlgorithm::computeSecondSpringConstant(const Framework::Node* const firstNode, 
						     const Framework::Node* const secondNode) 
  {
//////////// FB: Test 2 different spring constants still under constraction
    /* CFAUTOTRACE;
  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
    const CFreal firstNodeValue  = nodalStates[firstNode->getLocalID()] [m_monitorVarID+1];
    const CFreal secondNodeValue = nodalStates[secondNode->getLocalID()][m_monitorVarID+1];
    return std::abs(secondNodeValue - firstNodeValue);
    }*/
  /*
  cf_assert(m_monitorPhysVarID < m_pdata.size());
  // physical data arrays are computed on-the-fly from given nodal states 
  m_state->copyData(nodalStates[firstNode->getLocalID()]);
  getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
  const CFreal firstNodeValue  = m_pdata[m_monitorPhysVarID+1];
  m_state->copyData(nodalStates[secondNode->getLocalID()]);
  getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
  const CFreal secondNodeValue  = m_pdata[m_monitorPhysVarID];
  return std::abs(secondNodeValue - firstNodeValue);*/
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
  const CFuint nbNodes = nodes.size();
  Framework::DataHandle <CFreal> nodeDistance = socket_nodeDistance.getDataHandle();

  Common::SafePtr<Framework::TopologicalRegionSet> cells =
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode){
    //RSI computations
    if(m_tolerance != 0.){
      const std::string nsp = this->getMethodData().getNamespace();
      if(m_erreurG > m_tolerance/100.||  SubSystemStatusStack::getActive()->getNbIter()<1000){
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
	      if(insideRegion(nodes[iNode])){
		if(nbElemTypes==3 && nbDims==2){
		  assembleinRegionNode2DTriag(nodes[iNode]);
		  //assembleInnerNode(nodes[iNode]);
		}
		if(nbElemTypes==4 && nbDims==2){
		  assembleinRegionNode2DQuads(nodes[iNode]);
		}
		if(nbElemTypes==4 && nbDims==3){
		  assembleinRegionNode3DTet(nodes[iNode]);
		}
	      }
	      else{
		assembleInnerNode(nodes[iNode]);
	      }
	    }
	  }
	}
      }
    }

    // No RSI computations
    if(m_tolerance == 0.){
      if (!nodes[iNode]->isParUpdatable()){ 
	//do nothing
      }
      else{
	if(isNodeMovingInBoundary(nodes[iNode])  && insideRegion(nodes[iNode])==false){

	  assembleLockedNode(nodes[iNode]);
	}	
	else{ 
	  if(isNodeLocked(nodes[iNode])){
	    assembleLockedNode(nodes[iNode]);
	  }
	  else{
	    if(insideRegion(nodes[iNode]) ){
	      if(nbElemTypes==3 && nbDims==2){
		assembleinRegionNode2DTriag(nodes[iNode]);
	      }
	      if(nbElemTypes==4 && nbDims==2){
		assembleinRegionNode2DQuads(nodes[iNode]);	

	      }
	      if(nbElemTypes==4 && nbDims==3){
		assembleinRegionNode3DTet(nodes[iNode]);
	      }
	      if(nbElemTypes==8 && nbDims==3){
		assembleinRegionNode3DHexa(nodes[iNode]);
	      }
	    }
	    else{
	      assembleInnerNode(nodes[iNode]);
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
  CFLog(VERBOSE, "MeshFittingAlgorithm::Finishing assembleLinearSystem()\n");

}
/////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithm::computeState(const Framework::Node*  Node){
  CFAUTOTRACE;
  Framework::DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
    const CFreal firstNodeValue  = nodalStates[Node->getLocalID()] [m_monitorVarID];
    return (firstNodeValue);
  }
  
  cf_assert(m_monitorPhysVarID < m_pdata.size());
  // physical data arrays are computed on-the-fly from given nodal states 
  m_state->copyData(nodalStates[Node->getLocalID()]);
  getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
  const CFreal firstNodeValue  = m_pdata[m_monitorPhysVarID];
  return firstNodeValue;
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
      

//Joachim
//////////////////////////////////////////////////////////////////////////////

void MeshFittingAlgorithm::assembleMovingInBoundaryNode3DHexa(const Framework::Node* node){
   
  bool blocked = false;
  Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  const RealVector& nodeNormal = (m_mapNodeIDNormal[node->getLocalID()]); // get the normal of each moving in boundary node
  
  std::vector<Framework::Node*>::const_iterator it;
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
  std::vector<const Framework::Node*> sharedNodes;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    Framework::Node* neighborNode = *it;
    RealVector dist(nbDims);
    const bool neighborIsBoundary = m_boundaryNodes.find(neighborNode) != m_boundaryNodes.end();
    if (neighborIsBoundary){
      for (CFuint iDim=0; iDim<nbDims; ++iDim){
	dist[iDim] = (*neighborNode)[XX+iDim]-(*node)[XX+iDim];
      }
      if (dist.norm2() < m_equilibriumSpringLength*1 && blocked == false){
	blocked = true;
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
    
    CFreal sumOffDiagonalValues = 0.;
    const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node); // get neighbours
    std::vector<Framework::Node*>::const_iterator it;
    for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
      const Framework::Node* neighborNode = *it;
      const CFreal springConstant = computeSpringConstant(node, neighborNode); // spring constant
      
      RealVector springDirection(nbDims); // vector from the node to the neighbour
      for (CFuint iDim=0; iDim<nbDims; ++iDim){
	springDirection[iDim] = (*neighborNode)[XX+iDim] - (*node)[XX+iDim] ;
      }
      springDirection.normalize(); // normalized vector from the node to its neighbour
      RealVector projectedSpringDirection = springDirection - 
	MathTools::MathFunctions::innerProd(springDirection, nodeNormal)*nodeNormal;
      
      const CFreal dotProduct = MathTools::MathFunctions::innerProd(projectedSpringDirection, springDirection);
      
      const CFreal physicalSpringConstant =  truncateSpringConstant(springConstant)*std::abs(dotProduct);
      CFreal normalizedSpringConstant =  physicalSpringConstant;
      
      const CFreal stiff= normalizedSpringConstant;
      stiffness[node->getLocalID()]=stiff;
      sumOffDiagonalValues += stiff;
      
      const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
      for(CFuint iFreeDim=0; iFreeDim<freeDims.size(); ++iFreeDim){
	jacobMatrix->addValue(rowGlobalID+freeDims[iFreeDim], colGlobalID+freeDims[iFreeDim], stiff);
      }
    }
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

 //////////////////////////////////////////////////////////////////////////////
      
bool MeshFittingAlgorithm::isNodeLocked( Framework::Node* node)
{
  const bool isBoundary = m_boundaryNodes.find(node) != m_boundaryNodes.end();
  return ((isBoundary)) ;
}
     
 //////////////////////////////////////////////////////////////////////////////
      

bool MeshFittingAlgorithm::insideRegion(Framework::Node* node)
{
  bool inRegion=false;
  Framework::Node& currNode = *node;
  const CFuint nodeID = currNode.getLocalID();
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
 // FB : uncomment isNodeLocked(node)==false f you are not using nodal movememt interpolation
  if  ((nodeisAD[nodeID] == true)   && (isNodeLocked(node)==false)){
    inRegion=true;	    
  }
  return (inRegion);
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
RealVector  MeshFittingAlgorithm::computeIntersection(const  Framework::Node* const  a,const  Framework::Node* const  c,
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
///////////////////////////////

RealVector  MeshFittingAlgorithm::computeIntersectionQuad3D(const  Framework::Node* const  a,const  Framework::Node* const  c,
						      const  Framework::Node* const  b, const  Framework::Node* const  d){
  RealVector xyi(3); xyi =0.;
  xyi[0] = ((*a)[0]+(*b)[0]+(*c)[0]+(*d)[0])/4;
  xyi[1] = ((*a)[1]+(*b)[1]+(*c)[1]+(*d)[1])/4;
  xyi[2] = ((*a)[2]+(*b)[2]+(*c)[2]+(*d)[2])/4;
  return xyi;
}


 //////////////////////////////////////////////////////////////////////////////
  
CFreal MeshFittingAlgorithm::computeElementArea2dQuads(const  Framework::Node* const  a,const  Framework::Node* const  c,
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
//////////////////////////////////////////////////////////////////////////////
   
CFreal MeshFittingAlgorithm::computeConstantquads(const  RealVector xyi , 
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

  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    vector1[iDim]= (xyi)[iDim]-(*firstNode)[XX+iDim];
    vector2[iDim]= (xyi)[iDim]-(*secondNode)[XX+iDim];
    
    vector3[iDim]= (*thirdNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector4[iDim]= (*thirdNode)[XX+iDim]-(*secondNode)[XX+iDim];

    vector5[iDim]= (*fourthNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector6[iDim]= (*fourthNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }

  MathTools::MathFunctions::crossProd(vector1, vector2,res);
  const CFreal  area = (res.norm2()/2);

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

/////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithm::computeConstantquads3D(const  Framework::Node* const  firstNode,
						  const  Framework::Node* const  secondNode,
						  const Framework::Node* const  thirdNode,
						  const Framework::Node* const  fourthNode){
  
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
                 
  	
  
  RealVector vector3(nbDims); vector3 = 0.;
  RealVector vector4(nbDims); vector4 = 0.;
  RealVector vector5(nbDims); vector5 = 0.;
  RealVector vector6(nbDims); vector6 = 0.;
  CFreal torsionConstant;
  RealVector res(nbDims);
  RealVector res1(nbDims);
  RealVector res2(nbDims);

  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    
    vector3[iDim]= (*thirdNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector4[iDim]= (*thirdNode)[XX+iDim]-(*secondNode)[XX+iDim];

    vector5[iDim]= (*fourthNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector6[iDim]= (*fourthNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }

  MathTools::MathFunctions::crossProd(vector3, vector4,res);
  const CFreal  area1 = (res.norm2()/2);

  MathTools::MathFunctions::crossProd(vector5, vector6,res);
  const CFreal  area2 = (res.norm2()/2);



  CFreal k1 = (vector3.norm2()*vector3.norm2()*vector4.norm2()*vector4.norm2())/(4*area1*area1);

  CFreal k2 = (vector5.norm2()*vector5.norm2()*vector6.norm2()*vector6.norm2())/(4*area2*area2);
  // if( m_thetaMid ){
    // choose this option for a semi torsional spring analogy based on the middle angle only
  //torsionConstant = k;
  // }
  //else{
    // choose this option to stiffer the mesh or on a pave mesh
    torsionConstant = (std::max(k1,k2));
    // }
  return torsionConstant;
}

///////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithm::computeConstantquads3DElemArea(const  RealVector xyi , 
						  const  Framework::Node* const  firstNode,
						  const  Framework::Node* const  secondNode,
						  const Framework::Node* const  thirdNode,
						  const Framework::Node* const  fourthNode,
						  CFreal elementArea){
  
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
                 
  	
  
  RealVector vector1(nbDims); vector1 = 0.;
  RealVector vector2(nbDims); vector2 = 0.;
  RealVector vector3(nbDims); vector3 = 0.;
  RealVector vector4(nbDims); vector4 = 0.;
  RealVector vector5(nbDims); vector5 = 0.;
  RealVector vector6(nbDims); vector6 = 0.;
  CFreal torsionConstant;
  RealVector res(nbDims);
  RealVector res1(nbDims);
  RealVector res2(nbDims);

  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    vector1[iDim]= (xyi)[iDim]-(*firstNode)[XX+iDim];
    vector2[iDim]= (xyi)[iDim]-(*secondNode)[XX+iDim];
    
    vector3[iDim]= (*thirdNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector4[iDim]= (*thirdNode)[XX+iDim]-(*secondNode)[XX+iDim];

    vector5[iDim]= (*fourthNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector6[iDim]= (*fourthNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }
  MathTools::MathFunctions::crossProd(vector1, vector2,res);
  const CFreal  area = (res.norm2()/2);

  MathTools::MathFunctions::crossProd(vector3, vector4,res);
  const CFreal  area1 = (res.norm2()/2);

  MathTools::MathFunctions::crossProd(vector5, vector6,res);
  const CFreal  area2 = (res.norm2()/2);


  CFreal k = (area/elementArea)*(vector1.norm2()*vector2.norm2()*vector1.norm2()*vector2.norm2())/(4*area*area);

  CFreal k1 = (area1/elementArea)*(vector3.norm2()*vector3.norm2()*vector4.norm2()*vector4.norm2())/(4*area1*area1);

  CFreal k2 =(area2/elementArea)*(vector5.norm2()*vector5.norm2()*vector6.norm2()*vector6.norm2())/(4*area2*area2);

  if (k > 1000*k1 && k>1000*k2){
	cout << "k is the best" << endl;
  }

  // if( m_thetaMid ){
    // choose this option for a semi torsional spring analogy based on the middle angle only
  //torsionConstant = k;
  // }
  //else{
    // choose this option to stiffer the mesh or on a pave mesh
    torsionConstant =(std::max(k1,k2));
    // }
  return torsionConstant;
}

/////////////////////////////////////////////////////////////////////////////
CFreal MeshFittingAlgorithm::computetorsionConstant(const  Framework::Node* const  centralNode, 
						    const  Framework::Node* const  firstNode,
						    const  Framework::Node* const  secondNode){
  //----------
  //This part is using the states as the linear spring and is divided by the area to avoid zero volumes 

  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector vector1(3); vector1 = 0.;
  RealVector vector2(3); vector2 = 0.;
  RealVector res(3);
  CFreal torsionConstant;
  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    vector1[iDim]= (*centralNode)[XX+iDim]-(*firstNode)[XX+iDim];
    vector2[iDim]= (*centralNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }
  MathTools::MathFunctions::crossProd(vector1, vector2,res);
  const CFreal  area = (res.norm2()/2);
  CFreal k = (vector1.norm2()*vector2.norm2()*vector1.norm2()*vector2.norm2())/(4*area*area);
  torsionConstant = k;
  return torsionConstant;
}
//////////////////////////////////////////////////////////////////////////////
      
CFreal MeshFittingAlgorithm::computetorsionConstant3DSemiOnly(CFuint  nodeID1,CFuint nodeID2 ,
						    const  Framework::Node* const  firstNode,
						    const  Framework::Node* const  secondNode){

  //// semi torsional spring analogy for 3D tethrahedral

  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector projectedpoint(3); projectedpoint = 0.;
  RealVector node1(3); node1 = 0.;
  RealVector node2(3); node2 = 0.;
  RealVector edgeNode1Node2(3); edgeNode1Node2 = 0.;
  RealVector edgeNode2Node1(3); edgeNode2Node1 = 0.;
  RealVector Node2proj(3); Node2proj=0.;
  RealVector Node1proj(3); Node1proj=0.;
  RealVector edgenode2firstnode(3);edgenode2firstnode = 0.;
  RealVector edgenode1firstnode(3); edgenode1firstnode = 0.;
  RealVector projectedPointSecondNode(3) ; projectedPointSecondNode = 0.;
  RealVector projectedPointFirstNode(3) ; projectedPointFirstNode = 0.;
  RealVector res(3); res=0.; 

  RealVector vNeighborNode1(3) ; vNeighborNode1 = 0.;
  RealVector vNeighborNode2(3) ; vNeighborNode2 = 0.;
  RealVector vfirstNodeNeighbor(3) ;vfirstNodeNeighbor = 0.;
  RealVector normal(3);normal = 0. ; 
  RealVector length(3); length = 0.; 
  CFreal torsionConstant=0.;
  const std::vector<Framework::Node*>& neighboringNodesTest = m_edgeGraphN.getNeighborNodesOfNode(firstNode);
  std::vector<Framework::Node*>::const_iterator it;
  for(it=neighboringNodesTest.begin(); it != neighboringNodesTest.end(); ++it){
    const Framework::Node* neighborNodeTest = *it;
    if (neighborNodeTest->getLocalID() == nodeID1){
      for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
	node1[iDim] = (*neighborNodeTest)[XX+iDim];
      }
    }
    if (neighborNodeTest->getLocalID() == nodeID2){
      for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
	node2[iDim] = (*neighborNodeTest)[XX+iDim];
      }
    }
  }

    
  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    edgeNode1Node2[iDim]= node2[iDim]-node1[iDim];
    edgeNode2Node1[iDim]= node1[iDim]-node2[iDim];
    edgenode1firstnode[iDim]= node1[iDim]-(*firstNode)[XX+iDim];
    edgenode2firstnode[iDim]= node2[iDim]-(*firstNode)[XX+iDim];
    vNeighborNode1[iDim]=(*secondNode)[XX+iDim]-node1[iDim];
    vNeighborNode2[iDim]=(*secondNode)[XX+iDim]-node2[iDim];
    vfirstNodeNeighbor[iDim]=(*firstNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }

  edgeNode1Node2.normalize();
  edgeNode2Node1.normalize();

  if (edgenode1firstnode.norm2()>= edgenode2firstnode.norm2()){
    for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
      Node1proj[iDim]= MathTools::MathFunctions::innerProd(edgenode1firstnode,edgeNode1Node2)*edgeNode1Node2[iDim];
      projectedPointFirstNode[iDim]= Node1proj[iDim]+edgenode1firstnode[iDim];
      projectedPointSecondNode[iDim]=Node1proj[iDim]-vNeighborNode1[iDim];
    }
  }
  if (edgenode1firstnode.norm2()< edgenode2firstnode.norm2()){
    for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
      Node2proj[iDim]= MathTools::MathFunctions::innerProd(edgenode2firstnode,edgeNode1Node2)*edgeNode2Node1[iDim];
      projectedPointFirstNode[iDim]= Node2proj[iDim]+edgenode2firstnode[iDim];
      projectedPointSecondNode[iDim]=Node2proj[iDim]-vNeighborNode2[iDim];
    }
  }
  MathTools::MathFunctions::crossProd(projectedPointSecondNode,projectedPointFirstNode ,res);
  CFreal  area = (res.norm2()/2);
  const CFreal k =(projectedPointSecondNode.norm2()*projectedPointFirstNode.norm2()*(projectedPointSecondNode.norm2()*projectedPointFirstNode.norm2()))/(4*area*area);
  torsionConstant=(k);
  return torsionConstant;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MeshFittingAlgorithm::computetorsionConstant3D(CFuint  nodeID1,CFuint nodeID2 ,
						    const  Framework::Node* const  firstNode,
						    const  Framework::Node* const  secondNode){
  //// semi torsional spring analogy for 3D tethrahedral+OST

  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector projectedpoint(3); projectedpoint = 0.;
  RealVector node1(3); node1 = 0.;
  RealVector node2(3); node2 = 0.;
  RealVector edgeNode1Node2(3); edgeNode1Node2 = 0.;
  RealVector edgeNode2Node1(3); edgeNode2Node1 = 0.;

  RealVector edgenode2firstnode(3);edgenode2firstnode = 0.;
  RealVector edgenode1firstnode(3); edgenode1firstnode = 0.;
  RealVector projectedPointSecondNode(3) ; projectedPointSecondNode = 0.;
  RealVector projectedPointFirstNode(3) ; projectedPointFirstNode = 0.;
  RealVector res(3); res=0.; 
  RealVector res1(3); res1=0.;
  RealVector res2(3); res2=0.;
  RealVector Node2proj(3); Node2proj=0.;
  RealVector Node1proj(3); Node1proj=0.;

  RealVector vNeighborNode1(3) ; vNeighborNode1 = 0.;
  RealVector vNeighborNode2(3) ; vNeighborNode2 = 0.;
  RealVector vfirstNodeNeighbor(3) ;vfirstNodeNeighbor = 0.;
  RealVector normal(3);normal = 0. ; 
  //CFreal  firstNode_P = 0.;
  RealVector length(3); length = 0.; 
  //CFreal torsionConstant=0.;
  CFreal max=0.;
  const std::vector<Framework::Node*>& neighboringNodesTest = m_edgeGraphN.getNeighborNodesOfNode(firstNode);
  std::vector<Framework::Node*>::const_iterator it;
  for(it=neighboringNodesTest.begin(); it != neighboringNodesTest.end(); ++it){
    const Framework::Node* neighborNodeTest = *it;
    RealVector st(nbDims); st=0.;
    for(CFuint iDim=0; iDim<nbDims ; ++iDim){
      st[iDim] =(*firstNode)[XX+iDim]-(*neighborNodeTest)[XX+iDim];
    }
    if(max < st.norm2()){
      max=st.norm2(); 
    }
    if (neighborNodeTest->getLocalID() == nodeID1){
      for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
	node1[iDim] = (*neighborNodeTest)[XX+iDim];
      }
    }
    if (neighborNodeTest->getLocalID() == nodeID2){
      for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
	node2[iDim] = (*neighborNodeTest)[XX+iDim];
      }
    }
  }

    
  for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
    edgeNode1Node2[iDim]= node2[iDim]-node1[iDim];
    edgeNode2Node1[iDim]= node1[iDim]-node2[iDim];
    edgenode1firstnode[iDim]= node1[iDim]-(*firstNode)[XX+iDim];
    edgenode2firstnode[iDim]= node2[iDim]-(*firstNode)[XX+iDim];
    vNeighborNode1[iDim]=(*secondNode)[XX+iDim]-node1[iDim];
    vNeighborNode2[iDim]=(*secondNode)[XX+iDim]-node2[iDim];
    vfirstNodeNeighbor[iDim]=(*firstNode)[XX+iDim]-(*secondNode)[XX+iDim];
  }

  edgeNode1Node2.normalize();
  edgeNode2Node1.normalize();

  if (edgenode1firstnode.norm2()>= edgenode2firstnode.norm2()){
    for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
      Node1proj[iDim]= MathTools::MathFunctions::innerProd(edgenode1firstnode,edgeNode1Node2)*edgeNode1Node2[iDim];
      projectedPointFirstNode[iDim]= Node1proj[iDim]+edgenode1firstnode[iDim];
      projectedPointSecondNode[iDim]=Node1proj[iDim]-vNeighborNode1[iDim];
    }
  }
  if (edgenode1firstnode.norm2()< edgenode2firstnode.norm2()){
    for (CFuint iDim=0 ; iDim<nbDims ; iDim++){
      Node2proj[iDim]= MathTools::MathFunctions::innerProd(edgenode2firstnode,edgeNode1Node2)*edgeNode2Node1[iDim];
      projectedPointFirstNode[iDim]= Node2proj[iDim]+edgenode2firstnode[iDim];
      projectedPointSecondNode[iDim]=Node2proj[iDim]-vNeighborNode2[iDim];
    }
  }
  RealVector vNode1Neighbor(nbDims); vNode1Neighbor=0.;
  RealVector vNode2Neighbor(nbDims); vNode2Neighbor=0.;
  RealVector vNeighborFirst(nbDims); vNeighborFirst=0.;
  for(CFuint iDim =0; iDim<nbDims; ++iDim){
    vNode1Neighbor[iDim]=-vNeighborNode1[iDim];
    vNode2Neighbor[iDim]=-vNeighborNode2[iDim];
    vNeighborFirst[iDim]=-vfirstNodeNeighbor[iDim];
  }
 
  MathTools::MathFunctions::crossProd(projectedPointSecondNode,projectedPointFirstNode ,res);

  CFreal  area = (res.norm2()/2);
  RealVector d2(nbDims); d2=0.;
  RealVector d22(nbDims); d22=0.;
  RealVector n(nbDims); n=0.;


  for( CFuint iDim = 0; iDim< nbDims; ++iDim){
    d2[iDim] =(node1[iDim]-(*secondNode)[XX+iDim]);
    d22[iDim]=(node2[iDim]-(*secondNode)[iDim]);
  }

  MathTools::MathFunctions::crossProd(d2,d22 ,n);
  CFreal disToPlane = std::abs(MathTools::MathFunctions::innerProd(n,vNeighborFirst))/n.norm2();
  CFreal totalLenght = vfirstNodeNeighbor.norm2() + edgenode1firstnode.norm2()+ edgenode2firstnode.norm2();
  CFreal lambda1 = vfirstNodeNeighbor.norm2()/totalLenght;
  const CFreal kost1 = (1./disToPlane)/pow(lambda1,1);
  const CFreal d_nn_i =std::sqrt( vfirstNodeNeighbor.norm2()*vfirstNodeNeighbor.norm2()-disToPlane*disToPlane);
  const CFreal lambda2 = lambda1;
  const CFreal kost2 = (1./ d_nn_i)/pow(lambda2,1);


  const CFreal kost =(kost1+kost2)/2;

  const CFreal k =(projectedPointSecondNode.norm2()*projectedPointFirstNode.norm2()*(projectedPointSecondNode.norm2()*projectedPointFirstNode.norm2()))/(4*area*area);
  const  CFreal  torsionConstant=(kost+k);

  return torsionConstant;
}

 //////////////////////////////////////////////////////////////////////////////     
    
CFreal  MeshFittingAlgorithm::ComputeTvolume(CFuint  n1, 
					     CFuint  n2,
					     CFuint  n3,
					     CFuint  n4){
  
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  
  Framework::Node* node1;
  Framework::Node* node2;
  Framework::Node* node3;
  Framework::Node* node4;
  for(CFuint iNode=0; iNode<nodes.size(); ++iNode){
    if(nodes[iNode]->getLocalID()==n1){
      node1=nodes[iNode];
    }
    if(nodes[iNode]->getLocalID()==n2){
      node2=nodes[iNode];
    }
    if(nodes[iNode]->getLocalID()==n3){
      node3=nodes[iNode];
    }
    if(nodes[iNode]->getLocalID()==n4){
      node4=nodes[iNode];
    }
  }

  RealVector a(nbDims); a=0.;
  RealVector b(nbDims); b=0.;
  RealVector c(nbDims); c=0.;

  RealVector res(nbDims); res=0.;
  double dot =0.;
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
     a[iDim] =(*node1)[XX+iDim]-(*node2)[XX+iDim];
     b[iDim] =(*node1)[XX+iDim]-(*node3)[iDim];
     c[iDim] =(*node1)[XX+iDim]-(*node4)[iDim];

  }

  MathTools::MathFunctions::crossProd(a,b ,res);

  dot = MathTools::MathFunctions::innerProd(c, res);
  
  CFreal  volume = std::abs((1./6.)*dot);

  return volume;
}
      
//////////////////////////////////////////////////////////////////////////////
CFreal  MeshFittingAlgorithm::ComputeTFacesurface(  Framework::Node* node1, 
						    Framework::Node* node2,
						   Framework::Node* node3){
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  RealVector a(nbDims); a=0.;
  RealVector b(nbDims); b=0.;
  RealVector c(nbDims); c=0.;
  RealVector res(nbDims); res=0.;
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
     a[iDim] =(*node1)[XX+iDim]-(*node2)[XX+iDim];
     b[iDim] =(*node1)[XX+iDim]-(*node3)[iDim];
  }

  MathTools::MathFunctions::crossProd(a,b ,res);  
  CFreal  surface = res.norm2()/2.;


  return surface;
}
       
//////////////////////////////////////////////////////////////////////////////
 void MeshFittingAlgorithm::assembleinRegionNode2DTriag(const  Framework::Node* node){
 
   ////2D triangular 

  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  CFreal sumOffDiagonalValues = 0.;
  CFreal sum = 0.;
  std::vector<const Framework::Node*> sharedNodes;
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);

  std::vector<Framework::Node*>::const_iterator it;

  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    sharedNodes.clear();
    const Framework::Node* neighborNode = *it;
    const CFreal springConstant = computeSpringConstant(node,neighborNode);
    CFreal normalizedSpringConstant = truncateSpringConstant(springConstant);
    sumOffDiagonalValues +=normalizedSpringConstant ;
    CFreal kTorsion = 0.;
    const std::vector<Framework::Node*>& neighboringNodes1 = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    std::vector<Framework::Node*>::const_iterator it1;
    for(it1=neighboringNodes1.begin(); it1 != neighboringNodes1.end(); ++it1){
      Framework::Node* neighborNode1 = *it1;
      const std::vector<Framework::Node*>& neighboringNodesB = m_edgeGraph.getNeighborNodesOfNode(node);
      std::vector<Framework::Node*>::const_iterator it2;
      for(it2=neighboringNodesB.begin(); it2 != neighboringNodesB.end(); ++it2){
	Framework::Node* neighborNode2ndLoop = *it2;
	if (neighborNode1->getLocalID()==neighborNode2ndLoop->getLocalID()){
	  sharedNodes.push_back(neighborNode1);
	}
      }
    }
    for (CFuint b=0 ; b< sharedNodes.size() ; b++){
      std::vector<Framework::Node*>::const_iterator it3;
      for(it3=neighboringNodes.begin(); it3 != neighboringNodes.end(); ++it3){
	Framework::Node* neighborNode1B = *it3;
	if (neighborNode1B==sharedNodes[b]){
	  CFreal torsionConstant =computetorsionConstant(sharedNodes[b],neighborNode ,node); 
	  CFreal torsionTruncated  = (torsionConstant);
	  kTorsion += torsionTruncated;
	}
      }
    }
    sum +=kTorsion*normalizedSpringConstant;
    const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
    const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
    CFreal stiff =  normalizedSpringConstant+(normalizedSpringConstant*kTorsion);
    stiffness[node->getLocalID()]=stiff;
    for(CFuint iDim=0; iDim<nbDims; ++iDim){
      jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,stiff);
    }
    
  }
  //Diagonal value 
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue = -sum-sumOffDiagonalValues;    
    const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
  }
  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sum+sumOffDiagonalValues);
  }
  
 }
 //////////////////////////////////////////////////////////////////////////////
   
void MeshFittingAlgorithm::assembleinRegionNode2DQuads(const  Framework::Node* node){
  
      //2D quads
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
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
  std::vector<Framework::Node*>::const_iterator it;
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator itN1;
  std::vector<Framework::Node*>::const_iterator itSa;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    CFreal torsionConstant=0;
    std::vector<const Framework::Node*> sharedNodes;
    sharedNodes.clear();
    const  Framework::Node* neighborNode = *it;
    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
      const Framework::Node* neighborNodeOfN = *itN;
      for(itN1=neighboringNodes.begin(); itN1 != neighboringNodes.end(); ++itN1){
	const Framework::Node* neighborNodeNewloop = *itN1;
	if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
	  sharedNodes.push_back(neighborNodeNewloop);
	}
      }
    }
    if (sharedNodes.size()==4){
      for (mapIt it = ite.first; it != ite.second; ++it) {
	for(CFuint i=0; i<sharedNodes.size() ; ++i){
	  if (sharedNodes[i]->getLocalID() == it->second){
	    bool foundN = false;
	    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
	    cf_assert(foundN);
	    bool foundS = false;
	    std::pair<mapItS,mapItS > iteS=m_mapNodeNode1.find(sharedNodes[i]->getLocalID(), foundS);
	    cf_assert(foundS);
	    for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
	      for (mapItS itSe = iteS.first; itSe != iteS.second; ++itSe) {
		if (itNe->second == itSe->second){
		  for(itSa=neighboringNodes.begin(); itSa != neighboringNodes.end(); ++itSa){
		    if((*itSa)->getLocalID() == itNe->second && (*itSa)->getLocalID()!= node->getLocalID()){
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
      const CFreal springConstant =computeSpringConstant(node,neighborNode);
      const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant); 
      CFreal f = 1.;
      if (m_smoothSpringNetwork) {
      	//CFreal f = (7.-2.)/(0.02-0.01)*nodeDistance[node->getLocalID()]+ 1.2 - (7.-2.)/(0.02-0.01)*0.01;  //0.09
	// FB :test case dependent
      	CFreal f = (5.-1.)/(0.01-m_acceptableDistance)*nodeDistance[node->getLocalID()]+ 1. - (5.-1.)/(0.01-m_acceptableDistance)*m_acceptableDistance;  // 5
      }
      const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
      const CFreal stiff=normalizedSpringConstant+(pow(torsionConstant,1.)*pow(normalizedSpringConstant,f)); 
      //const CFreal stifNormalized = (stiff);
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

///////////////////////////////////////////////////////////////////////      

  void MeshFittingAlgorithm::assembleinRegionNode2DQuadsGeoBased(const  Framework::Node* node){
  ////////// Try: 2D based on geometry only 
      //2D quads
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
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
  std::vector<Framework::Node*>::const_iterator it;
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator itN1;
  std::vector<Framework::Node*>::const_iterator itSa;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    CFreal torsionConstant=0;
    std::vector<const Framework::Node*> sharedNodes;
    sharedNodes.clear();
    const  Framework::Node* neighborNode = *it;
    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
      const Framework::Node* neighborNodeOfN = *itN;
      for(itN1=neighboringNodes.begin(); itN1 != neighboringNodes.end(); ++itN1){
	const Framework::Node* neighborNodeNewloop = *itN1;
	if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
	  sharedNodes.push_back(neighborNodeNewloop);
	}
      }
    }
    if (sharedNodes.size()==4){
      for (mapIt it = ite.first; it != ite.second; ++it) {
	for(CFuint i=0; i<sharedNodes.size() ; ++i){
	  if (sharedNodes[i]->getLocalID() == it->second){
	    bool foundN = false;
	    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
	    cf_assert(foundN);
	    bool foundS = false;
	    std::pair<mapItS,mapItS > iteS=m_mapNodeNode1.find(sharedNodes[i]->getLocalID(), foundS);
	    cf_assert(foundS);
	    for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
	      for (mapItS itSe = iteS.first; itSe != iteS.second; ++itSe) {
		if (itNe->second == itSe->second){
		  for(itSa=neighboringNodes.begin(); itSa != neighboringNodes.end(); ++itSa){
		    if((*itSa)->getLocalID() == itNe->second && (*itSa)->getLocalID()!= node->getLocalID()){
		      RealVector d_node_neighbor(nbDims); d_node_neighbor=0.;
		      RealVector d_neighbor_itSa(nbDims); d_neighbor_itSa=0.;
		      for (CFuint iDim=0; iDim<nbDims; ++iDim){
			d_node_neighbor[iDim] = (*node)[XX+iDim] - (*neighborNode)[XX+iDim];
			d_neighbor_itSa[iDim] = (*neighborNode)[XX+iDim] - (*(*itSa))[XX+iDim];
		      }
		      // Activate to put ST on the smallest edge and Linear spring on the longest edge : Recommemded for high aspect ratio meshes.
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
      sum +=torsionConstant;
      const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
      const CFreal stiff=torsionConstant;
      stiffness[node->getLocalID()]=stiff;
      for(CFuint iDim=0; iDim<nbDims; ++iDim){
	jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,stiff);
      }
      
    }
  }
  
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue = -sum;
    const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
  }
  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue2 = sum;
    const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(diagValue2);
  }
 }

//////////////

bool MeshFittingAlgorithm::areCoplanar(const  Framework::Node* const  firstNode,
				       const  Framework::Node* const  secondNode,
				       const Framework::Node* const  thirdNode,
				       const Framework::Node* const  fourthNode){
						  
   CFuint x1 = (*firstNode)[0]; CFuint x2 = (*secondNode)[0]; CFuint x3 = (*thirdNode)[0]; CFuint x = (*fourthNode)[0];
   CFuint y1 = (*firstNode)[1]; CFuint y2 = (*secondNode)[1]; CFuint y3 = (*thirdNode)[1]; CFuint y = (*fourthNode)[1];
   CFuint z1 = (*firstNode)[2]; CFuint z2 = (*secondNode)[2]; CFuint z3 = (*thirdNode)[2]; CFuint z = (*fourthNode)[2];
   int a1 = x2 - x1 ; 
   int b1 = y2 - y1 ; 
   int c1 = z2 - z1 ; 
   int a2 = x3 - x1 ; 
   int b2 = y3 - y1 ; 
   int c2 = z3 - z1 ; 
   int a = b1 * c2 - b2 * c1 ; 
   int b = a2 * c1 - a1 * c2 ; 
   int c = a1 * b2 - b1 * a2 ; 
   int d = (- a * x1 - b * y1 - c * z1) ; 
   bool areCoplanar = false;     
    
   if(a * x + b * y + c * z + d < 1e-8){
     areCoplanar = true;
   }
   else{
     //cout << "not copl" << endl;
    }
   return areCoplanar;
}
/////////////

void MeshFittingAlgorithm::assembleinRegionNode3DHexa(const  Framework::Node* node){
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
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
  std::vector<Framework::Node*>::const_iterator it;
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator itN1;
  std::vector<Framework::Node*>::const_iterator itSa;
  
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    RealVector d_node_neighbor(nbDims); d_node_neighbor=0.;
    CFreal torsionConstant=0;
    std::vector<const Framework::Node*> sharedNodes;
    sharedNodes.clear();
    const  Framework::Node* neighborNode = *it;
    bool foundN = false;
    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
    cf_assert(foundN);
    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
      const Framework::Node* neighborNodeOfN = *itN;
      for(itN1=neighboringNodes.begin(); itN1 != neighboringNodes.end(); ++itN1){
	const Framework::Node* neighborNodeNewloop = *itN1;
	if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
	  sharedNodes.push_back(neighborNodeNewloop);
	}
      }
    }
    if (sharedNodes.size() == 16){ //edge-connected nodes
      for (CFuint iDim=0; iDim<nbDims; ++iDim){
	d_node_neighbor[iDim] = (*node)[XX+iDim] - (*neighborNode)[XX+iDim];
      }
      for (mapIt it = ite.first; it != ite.second; ++it) {
	for(CFuint i=0; i<sharedNodes.size() ; ++i){
	  if (sharedNodes[i]->getLocalID() == it->second){
	    bool foundN = false;	
	    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
	    cf_assert(foundN);
	    
	    bool foundS = false;
	    std::pair<mapItS,mapItS > iteS=m_mapNodeNode1.find(sharedNodes[i]->getLocalID(), foundS);
	    cf_assert(foundS);
	    for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
	      for (mapItS itSe = iteS.first; itSe != iteS.second; ++itSe) {
		if (itNe->second == itSe->second){
		  for(itSa=neighboringNodes.begin(); itSa != neighboringNodes.end(); ++itSa){
		    if((*itSa)->getLocalID() == itNe->second && (*itSa)->getLocalID()!= node->getLocalID()){
		      CFreal k3 = computetorsionConstant(*itSa, node, neighborNode);
		      CFreal k4 = computetorsionConstant(sharedNodes[i], node, neighborNode);
		      //torsionConstant += std::max(k3,k4); // to uncomment if torsional constants on faces are needed
		    }
		  }
		}	
	      }
	    }
	  }
	}
      }
      bool found = false;
      std::pair<mapIt,mapIt > iteCell=m_mapNodeCell1.find(node->getLocalID(), found);
      cf_assert(found);
      bool foundN = false;
      std::pair<mapItN,mapItN > iteNCell=m_mapNodeCell1.find(neighborNode->getLocalID(), foundN);
      cf_assert(foundN);
      std::vector<CFuint> cellsInCommon;
      for (mapIt it = iteCell.first; it != iteCell.second; ++it){
	for (mapItN itn = iteNCell.first; itn != iteNCell.second; ++itn){
	  if (it->second == itn->second){
	    cellsInCommon.push_back(it->second);
	  }
	}
      }
      for (CFuint m = 0; m < cellsInCommon.size(); ++m){
	bool foundCell = false;
	std::pair<mapIt,mapIt > itCell =m_mapCellNode1.find(cellsInCommon[m], foundCell);
	cf_assert(foundCell);
	std::vector<CFuint> oppositeNodes;
	for (mapIt itcell = itCell.first; itcell!= itCell.second; ++itcell){
	  bool foundNode = false;
	  for (mapIt it = ite.first; it != ite.second; ++it){
	    if (it->second == itcell->second){
	      foundNode = true;
	    }
	  }
	  for (mapItN itn  = iteN.first; itn != iteN.second; ++itn){
	    if (itn->second == itcell->second){
	      foundNode = true;
	    }
	  }
	  if (foundNode == false){
	    oppositeNodes.push_back(itcell->second);
	  }
	}
	bool foundO = false;
	std::pair<mapItS,mapItS > iteO=m_mapNodeNode1.find(oppositeNodes[0], foundO);
	cf_assert(foundO);
	bool isNodeL = false;
	for (mapItS its = iteO.first; its!=iteO.second; ++its){
	  for (mapIt it = ite.first; it != ite.second; ++it){
	    if (it->second == its->second){
	      isNodeL = true;
	    }
	  }
	}
	CFuint nodeL; 
	CFuint nodeK;
	if (isNodeL){
	  nodeL = oppositeNodes[0];
	  nodeK = oppositeNodes[1];
	}
	else{
	  nodeL = oppositeNodes[1];
	  nodeK = oppositeNodes[0];
	}
	const  Framework::Node* nodeNextI;
	const  Framework::Node* nodeNextJ;
	for(itSa=neighboringNodes.begin(); itSa != neighboringNodes.end(); ++itSa){
	  if((*itSa)->getLocalID() == nodeL){
	    nodeNextI = *itSa;
	  }
	  if((*itSa)->getLocalID() == nodeK){
	    nodeNextJ = *itSa;
	  }
	}
	CFreal k1 = computetorsionConstant(nodeNextJ, node, neighborNode);
	CFreal k2 = computetorsionConstant(nodeNextI, node, neighborNode);
	torsionConstant += std::max(k1,k2);			
      }
      
      const CFreal springConstant =computeSpringConstant(node,neighborNode);
      const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
      CFreal f = 1;
      if (m_smoothSpringNetwork) {
	//CFreal f = (7.-2.)/(0.02-0.01)*nodeDistance[node->getLocalID()]+ 1.2 - (7.-2.)/(0.02-0.01)*0.01;  //0.09
	// FB :test case dependent
	//cout << "smooth network" << endl;
	f =  1 + ((5.-1.)/(0.04-m_acceptableDistance))*nodeDistance[node->getLocalID()] - ((5.-1.)/(0.04-m_acceptableDistance))*m_acceptableDistance;  // 5
      }
      const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
      const CFreal stiff=normalizedSpringConstant;//    +(pow(torsionConstant,1.)*pow(normalizedSpringConstant,f)); 
      sum += stiff;
      stiffness[node->getLocalID()]=stiff;
      for(CFuint iDim=0; iDim<nbDims; ++iDim){
	jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,stiff);
      }
    }
  }
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue = -sum;
    const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
  }
  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue2 = sum;
    const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(diagValue2);
  }
  CFLog(VERBOSE, "MeshFittingAlgorithm::updateNodePositionshexaFinished()\n"); 

}




  //////////////////////////////
 void MeshFittingAlgorithm::smooth(const  Framework::Node* node){
  CFAUTOTRACE;
      //2D quads
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
  Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
  //cout<< nodes.size() << endl;
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle();
  const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  CFreal sumOffDiagonalValues = 0.;
  typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
  typedef MapNodeNode::MapIterator mapIt;
  typedef MapNodeNode::MapIterator mapItN;
  typedef MapNodeNode::MapIterator mapItS;
  bool found = false;
  bool blocked = false;
  RealVector term2(nbDims) ; term2 = 0.;
  RealVector term1(nbDims) ; term1 = 0.;
  
  std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
  cf_assert(found);
  const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
  std::vector<Framework::Node*>::const_iterator it;
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator itN1;
  std::vector<Framework::Node*>::const_iterator itSa;
  CFreal stiff;
  RealVector NN_Stiff(nodes.size()); NN_Stiff = 0.;
  std::vector<CFuint> nn_ID(4);
  CFuint it_nn = 0;
  RealVector num(nbDims); num = 0.; 
  RealVector den(nbDims); den = 0.;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    std::vector<const Framework::Node*> sharedNodes;
    sharedNodes.clear();
    const  Framework::Node* neighborNode = *it;
    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
      const Framework::Node* neighborNodeOfN = *itN;
      for(itN1=neighboringNodes.begin(); itN1 != neighboringNodes.end(); ++itN1){
	const Framework::Node* neighborNodeNewloop = *itN1;
	if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
	  sharedNodes.push_back(neighborNodeNewloop);
	}
      }
    }
    if (sharedNodes.size()==4){
      for (mapIt it = ite.first; it != ite.second; ++it) {
	for(CFuint i=0; i<sharedNodes.size() ; ++i){
	  if (sharedNodes[i]->getLocalID() == it->second){
	    bool foundN = false;
	    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
	    cf_assert(foundN); 
	    RealVector dist(nbDims); dist = 0.;
	    RealVector X(nbDims); X[0]=1.; X[1]=0.;
	    RealVector Y(nbDims); Y[0]=0.; Y[1]=1.;
	    RealVector distance(nbDims); distance = 0.;
	    RealVector unitDirection(nbDims); unitDirection = 0.;

	    for (CFuint iDim=0; iDim<nbDims; ++iDim){
	      dist[iDim] = (*neighborNode)[XX+iDim]-(*node)[XX+iDim];
	    }
	    for (CFuint iDim=0; iDim<nbDims; ++iDim){
	      unitDirection[iDim] = dist[iDim]/dist.norm2();
	    }
	    distance[0] =dist[0]- 100.*(MathTools::MathFunctions::innerProd(unitDirection,X))*m_equilibriumSpringLength; 
	    distance[1] =dist[1]- 100.*(MathTools::MathFunctions::innerProd(unitDirection,Y))*m_equilibriumSpringLength; 
	    const CFreal springConstant =computeSpringConstant(node,neighborNode);
	    const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
	    
	    for (CFuint iDim=0; iDim<nbDims; ++iDim){
		CFreal d = dist.norm2()/(m_equilibriumSpringLength*1000.);
	        //num[iDim] = num[iDim] + (1.-pow(d,4.))*exp(-pow(d,4.))*unitDirection[iDim];
		num[iDim] = num[iDim]+ (*neighborNode)[XX+iDim]*(dist.norm2()-100.*m_equilibriumSpringLength); 
	        den[iDim] = den[iDim]+ (dist.norm2()-100.*m_equilibriumSpringLength);   
	    }
	  }
	}
      }
    }
  }
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
    if(node->getLocalID() == nodes[iNode]->getLocalID() && nodeDistance[nodes[iNode]->getLocalID()] < .02 ){ 
      Framework::Node& currNode = *nodes[iNode];
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
	CFreal f=(.5/(0.02-0.01)) * nodeDistance[nodes[iNode]->getLocalID()] - (0.5*0.01/(0.02-0.01)) + 0.5;
	currNode[XX+iDim] = num[iDim]/den[iDim]*0.5 + currNode[XX+iDim]*0.5 ;
      }
    }
  }
 }
      //////////////////////////////////////////////////////////////////////////////////

 void MeshFittingAlgorithm::assembleinRegionNode3DTet(const  Framework::Node* node){

   //////3D tethrahedral 
   Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
   Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
   Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
   Common::SafePtr<Framework::TopologicalRegionSet> cells = 
     Framework::MeshDataStack::getActive()->getTrs("InnerCells");
   const CFuint nbCells = cells->getLocalNbGeoEnts();
   Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
   const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
   const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
   const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
   CFreal sumOffDiagonalValues = 0.;
   CFreal sum = 0.;
   typedef CFMultiMap<CFuint, CFuint> MapCellNode;
   typedef MapCellNode::MapIterator mapItc;
   const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
   std::vector<Framework::Node*>::const_iterator it;
   std::vector<Framework::Node*>::const_iterator itS;
   for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
     CFuint comp=0;
     CFreal torsionConstant=0;
     RealVector node1(3); node1 = 0.;
     RealVector node2(3); node2 = 0.;
     const  Framework::Node* neighborNode = *it;
     for (CFuint iCell=0; iCell<nbCells ; ++iCell){
       std::vector<CFuint> sharedNode;
       bool isThere = false;
       bool isThere2 = false;
       sharedNode.clear();
       bool foundC = false;
       std::pair<mapItc,mapItc > iteC=m_mapCellNode1.find(iCell, foundC);
       cf_assert(foundC);
       for (mapItc it2 = iteC.first; it2 != iteC.second; ++it2){
	 if (node->getLocalID()==it2->second){
	   isThere = true;
	 }
	 if (neighborNode->getLocalID()==it2->second){
	   isThere2 = true;
	 }
       }
       if (isThere==true && isThere2==true){
	 for (mapItc it2 = iteC.first; it2 != iteC.second; ++it2){
	   if((it2->second != neighborNode->getLocalID()) && (it2->second != node->getLocalID())){
	     sharedNode.push_back(it2->second);
	   }
	 }
       }
       
       if( sharedNode.size()==2){     
	torsionConstant +=computetorsionConstant3D(sharedNode[0],sharedNode[1],node,neighborNode);
       }
     } 
     const CFreal springConstant = computeSpringConstant(node,neighborNode);
     const CFreal normalizedSpringConstant = truncateSpringConstant(springConstant);
     sumOffDiagonalValues +=normalizedSpringConstant ;
     sum +=(torsionConstant)*normalizedSpringConstant;
     const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
     const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
     const CFreal stif =normalizedSpringConstant*torsionConstant+normalizedSpringConstant;
     stiffness[node->getLocalID()]=stif;
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,stif);
     }  
   }
   for(CFuint iDim=0; iDim<nbDims; ++iDim){
     CFreal diagValue = -sumOffDiagonalValues-sum;
     const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
     jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
   }
   //Right hand side
   for(CFuint iDim=0; iDim<nbDims; ++iDim){
    const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues+sum);
   }
 }


      
      ////////////////////////////// 
void MeshFittingAlgorithm::saveOldMeshProperties(){
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();

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
}
      /////////////////////////:
void MeshFittingAlgorithm::interpolateNewMeshProperties(){
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
    


      ///////////////////////////////:

void MeshFittingAlgorithm::stateInterpolator(){
  CFAUTOTRACE;
      //2D quads
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle< std::vector< Framework::State*> >  stencil = socket_stencil.getDataHandle();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  for(CFuint iState = 0; iState < states.size(); ++iState) {
    RealVector interpolatedStates(totalNbEqs); interpolatedStates=0.;
    RealVector oldState(totalNbEqs); oldState=0.;

    CFreal invDistance = 0.;
    const RealVector& stateCoord = states[iState]->getCoordinates();
    const CFuint stencilSize = stencil[iState].size();
    //cout<<" stencilSize " <<stencilSize << endl;
    // loop over the neighbor cells belonging to the chosen stencil
    for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
      oldState[iVar] = (*states[iState])[iVar];
    }
    for(CFuint in = 0; in < stencilSize; ++in) {
      const RealVector& stencilCoord = stencil[iState][in]->getCoordinates();
      //cout << "stencilCoord " << stencilCoord << endl;
      CFreal  stateTostencil = MathTools::MathFunctions::getDistance(stencilCoord,stateCoord );
      //cout << " stateTostencil "<<stateTostencil << endl;
      ///WTF needs to be fixed 
      //if(!oldState->isGhost()) {
      invDistance += 1./stateTostencil;

      for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
	interpolatedStates[iVar] += ((*stencil[iState][in])[iVar])*invDistance ;
	//cout << " ((*stencil[iState][in])[iVar]) " <<  ((*stencil[iState][in])[iVar])  <<endl;
	//cout << " invDis " << invDistance<<endl; ;
      }
    } 
    //}
    for(CFuint iVar = 0; iVar <totalNbEqs ; ++iVar) {
      //cout<< "Ivar " << iVar << endl;
      // cout<< "Oldstates[iState][iVar] " << (*states[iState])[iVar]<< endl;
      (*states[iState])[iVar] = interpolatedStates[iVar]/invDistance;
      //cout<< "Newstates[iState][iVar] " << (*states[iState])[iVar]<< endl;
    }
  } 
}
      








      
//////////////////////////////////////////////////////////////////////////////
 
 void MeshFittingAlgorithm::assembleInnerNode(const Framework::Node* node){
   // quads to be handeled
   
    SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
     getTrs("InnerCells");
   const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
   Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
   Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
   Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
   Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
   const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
   const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
   const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
   ////2D triangular or 3D tetrahedral
   if( (nbElemTypes==3 && nbDims==2) || (nbElemTypes==4 && nbDims==3) ){
     CFreal sumOffDiagonalValues = 0.;
     std::vector<const Framework::Node*> sharedNodes;
     const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
     std::vector<Framework::Node*>::const_iterator it;
     for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
       const Framework::Node* neighborNode = *it;
       const CFreal springConstant = computeSpringConstant(node,neighborNode);
       RealVector d_node_neighbor(nbDims); d_node_neighbor=0.;
       for (CFuint iDim=0; iDim<nbDims; ++iDim){
	 d_node_neighbor[iDim] = (*node)[XX+iDim] - (*neighborNode)[XX+iDim];
       }
       CFreal normalizedSpringConstant =1./d_node_neighbor.norm2()*truncateSpringConstant(springConstant);
       sumOffDiagonalValues += normalizedSpringConstant;
       const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
       const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
       CFreal stiff = normalizedSpringConstant;
       stiffness[node->getLocalID()]=stiff;
       for(CFuint iDim=0; iDim<nbDims; ++iDim){
	 jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,normalizedSpringConstant);
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
       rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues);
     }
   }
   
   //2D quads
   if( nbElemTypes==4 && nbDims==2){
     CFreal sumOffDiagonalValues = 0.;
     typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
     typedef MapNodeNode::MapIterator mapIt;
     typedef MapNodeNode::MapIterator mapItN;
     typedef MapNodeNode::MapIterator mapItS;
     bool found = false;
     std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
     cf_assert(found);
     const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
     std::vector<Framework::Node*>::const_iterator it;
     std::vector<Framework::Node*>::const_iterator itN;
     std::vector<Framework::Node*>::const_iterator itN1;
     std::vector<Framework::Node*>::const_iterator itSa;
     
     for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
       const  Framework::Node* neighborNode = *it;
       for (mapIt it = ite.first; it != ite.second; ++it) {
	 if (neighborNode->getLocalID() == it->second){
	   const CFreal springConstant =computeSpringConstant(node,neighborNode);
	   const  CFreal normalizedSpringConstant = truncateSpringConstant(springConstant);
	   sumOffDiagonalValues +=normalizedSpringConstant ;
	    stiffness[node->getLocalID()]= normalizedSpringConstant;

	   const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
	   const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
	   for(CFuint iDim=0; iDim<nbDims; ++iDim){
	     jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,normalizedSpringConstant);
	   }
	 }
    }
     }
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
       jacobMatrix->addValue(globalID+iDim, globalID+iDim, -sumOffDiagonalValues);
     }
     //Right hand side
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
       rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues);
     }
   
   }

   // 3D hexa
   if (nbElemTypes==8 && nbDims==3){
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
	  std::vector<Framework::Node*>::const_iterator it;
	  std::vector<Framework::Node*>::const_iterator itN;
	  std::vector<Framework::Node*>::const_iterator itN1;
	  std::vector<Framework::Node*>::const_iterator itSa;
	  
	  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
	    RealVector d_node_neighbor(nbDims); d_node_neighbor=0.;
	    CFreal torsionConstant=0;
	    std::vector<const Framework::Node*> sharedNodes;
	    sharedNodes.clear();
	    const  Framework::Node* neighborNode = *it;
	    bool foundN = false;
	    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
	    cf_assert(foundN);
	    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
	    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
	      const Framework::Node* neighborNodeOfN = *itN;
	      for(itN1=neighboringNodes.begin(); itN1 != neighboringNodes.end(); ++itN1){
		const Framework::Node* neighborNodeNewloop = *itN1;
		if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
		  sharedNodes.push_back(neighborNodeNewloop);
		}
	      }
	    }
	    if (sharedNodes.size() == 16){
		 const CFreal springConstant =computeSpringConstant(node,neighborNode);
		 const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
		 const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
		 const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
		 const CFreal stiff=normalizedSpringConstant;
		 sum += stiff;
		 stiffness[node->getLocalID()]=stiff;
		 for(CFuint iDim=0; iDim<nbDims; ++iDim){
		     jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,stiff);
		 }
	    }
	  }
	  for(CFuint iDim=0; iDim<nbDims; ++iDim){
	     CFreal diagValue = -sum;
	     const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
	     jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
	  }
	  //Right hand side
	  for(CFuint iDim=0; iDim<nbDims; ++iDim){
	     CFreal diagValue2 = sum;
	     const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
	     rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(diagValue2);
	  }
   }	

 }

   ///////////////////////////////// For test 

      //2D quads   For Qarman 
   /* Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
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
  std::vector<Framework::Node*>::const_iterator it;
  std::vector<Framework::Node*>::const_iterator itN;
  std::vector<Framework::Node*>::const_iterator itN1;
  std::vector<Framework::Node*>::const_iterator itSa;
  CFuint f=0;
  for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
    CFreal torsionConstant=0;
    std::vector<const Framework::Node*> sharedNodes;
    sharedNodes.clear();
    const  Framework::Node* neighborNode = *it;
    const std::vector<Framework::Node*>& neighboringNodesOfN = m_edgeGraphN.getNeighborNodesOfNode(neighborNode);
    for(itN=neighboringNodesOfN.begin(); itN != neighboringNodesOfN.end(); ++itN){
      const Framework::Node* neighborNodeOfN = *itN;
      for(itN1=neighboringNodes.begin(); itN1 != neighboringNodes.end(); ++itN1){
	const Framework::Node* neighborNodeNewloop = *itN1;
	if(neighborNodeOfN->getLocalID()==neighborNodeNewloop->getLocalID()){
	  sharedNodes.push_back(neighborNodeNewloop);
	}
      }
    }
    if (sharedNodes.size()==4){
      for (mapIt it = ite.first; it != ite.second; ++it) {
	for(CFuint i=0; i<sharedNodes.size() ; ++i){
	  if (sharedNodes[i]->getLocalID() == it->second){
	    bool foundN = false;
	    std::pair<mapItN,mapItN > iteN=m_mapNodeNode1.find(neighborNode->getLocalID(), foundN);
	    cf_assert(foundN);
	    bool foundS = false;
	    std::pair<mapItS,mapItS > iteS=m_mapNodeNode1.find(sharedNodes[i]->getLocalID(), foundS);
	    cf_assert(foundS);
	    for (mapItN itNe = iteN.first; itNe != iteN.second; ++itNe) {
	      for (mapItS itSe = iteS.first; itSe != iteS.second; ++itSe) {
		if (itNe->second == itSe->second){
		  for(itSa=neighboringNodes.begin(); itSa != neighboringNodes.end(); ++itSa){
		    if((*itSa)->getLocalID() == itNe->second && (*itSa)->getLocalID()!= node->getLocalID()){
		      // Activate to put ST on the smallest edge and Linear spring on the longest edge 
		      //RealVector d_node_neighbor(nbDims); d_node_neighbor=0.;
		      //RealVector d_neighbor_itSa(nbDims); d_neighbor_itSa=0.;
		      //for (CFuint iDim=0; iDim<nbDims; ++iDim){
			//d_node_neighbor[iDim] = (*node)[XX+iDim] - (*neighborNode)[XX+iDim];
			//d_neighbor_itSa[iDim] = (*neighborNode)[XX+iDim] - (*(*itSa))[XX+iDim];
			//}
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
      const CFreal springConstant =computeSpringConstant(node,neighborNode);
      const  CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
      sumOffDiagonalValues +=((normalizedSpringConstant));
      const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
      const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
      const CFreal stiff=((normalizedSpringConstant));
      stiffness[node->getLocalID()]=stiff;
      for(CFuint iDim=0; iDim<nbDims; ++iDim){
	jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,stiff);
      }
      f=f+1;
    }
  }
  
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue = -sumOffDiagonalValues;
    const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
    jacobMatrix->addValue(globalID+iDim, globalID+iDim, diagValue);
  }
  //Right hand side
  for(CFuint iDim=0; iDim<nbDims; ++iDim){
    CFreal diagValue2 = sumOffDiagonalValues;
    const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
    rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(diagValue2);
  }
  }*/


////////////// For test the second monitor variable
 
 void MeshFittingAlgorithm::assembleInnerNodeSecond(const Framework::Node* node){
   // quads to be handeled
   
    SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
     getTrs("InnerCells");
   const CFuint nbElemTypes = cells->getNbNodesInGeo(0);
   Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = socket_nodes.getDataHandle();
   Framework::DataHandle<CFreal> stiffness = socket_stiffness.getDataHandle();
   Framework::DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
   Common::SafePtr<Framework::LSSMatrix> jacobMatrix = m_lss->getMatrix();
   const Framework::LSSIdxMapping& idxMapping = m_lss->getLocalToGlobalMapping();
   const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
   const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
   ////2D triangular or 3D tetrahedral
   if( (nbElemTypes==3 && nbDims==2) || (nbElemTypes==4 && nbDims==3) ){
     CFreal sumOffDiagonalValues = 0.;
     std::vector<const Framework::Node*> sharedNodes;
     const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
     std::vector<Framework::Node*>::const_iterator it;
     for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
       const Framework::Node* neighborNode = *it;
       const CFreal springConstant = computeSecondSpringConstant(node,neighborNode);
       CFreal normalizedSpringConstant =truncateSpringConstant(springConstant);
       sumOffDiagonalValues += normalizedSpringConstant;
       const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
       const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
       CFreal stiff = normalizedSpringConstant;
       stiffness[node->getLocalID()]=stiff;
       for(CFuint iDim=0; iDim<nbDims; ++iDim){
	 jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,normalizedSpringConstant);
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
       rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues);
     }
   }
   
   //2D quads
   if( nbElemTypes==4 && nbDims==2){
     CFreal sumOffDiagonalValues = 0.;
     typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
     typedef MapNodeNode::MapIterator mapIt;
     typedef MapNodeNode::MapIterator mapItN;
     typedef MapNodeNode::MapIterator mapItS;
     bool found = false;
     std::pair<mapIt,mapIt > ite=m_mapNodeNode1.find(node->getLocalID(), found);
     cf_assert(found);
     const std::vector<Framework::Node*>& neighboringNodes = m_edgeGraph.getNeighborNodesOfNode(node);
     std::vector<Framework::Node*>::const_iterator it;
     std::vector<Framework::Node*>::const_iterator itN;
     std::vector<Framework::Node*>::const_iterator itN1;
     std::vector<Framework::Node*>::const_iterator itSa;
     
     for(it=neighboringNodes.begin(); it != neighboringNodes.end(); ++it){
       const  Framework::Node* neighborNode = *it;
       for (mapIt it = ite.first; it != ite.second; ++it) {
	 if (neighborNode->getLocalID() == it->second){
	   const CFreal springConstant =computeSecondSpringConstant(node,neighborNode);
	   const  CFreal normalizedSpringConstant = truncateSpringConstant(springConstant);
	   sumOffDiagonalValues +=normalizedSpringConstant ;
	    stiffness[node->getLocalID()]= normalizedSpringConstant;

	   const CFuint rowGlobalID = idxMapping.getRowID(node->getLocalID())*nbDims;
	   const CFuint colGlobalID = idxMapping.getColID(neighborNode->getLocalID())*nbDims;
	   for(CFuint iDim=0; iDim<nbDims; ++iDim){
	     jacobMatrix->addValue(rowGlobalID+iDim, colGlobalID+iDim,normalizedSpringConstant);
	   }
	 }
    }
     }
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       const CFuint globalID = idxMapping.getRowID(node->getLocalID())*nbDims;  
       jacobMatrix->addValue(globalID+iDim, globalID+iDim, -sumOffDiagonalValues);
     }
     //Right hand side
     for(CFuint iDim=0; iDim<nbDims; ++iDim){
       const CFreal equilibriumLength = m_equilibriumSpringLength*m_ratioBoundaryToInnerEquilibriumSpringLength;
       rhs[node->getLocalID()*totalNbEqs+XX+iDim] = equilibriumLength*(sumOffDiagonalValues);
     }
   
   }
 }
/////////////

void MeshFittingAlgorithm::computeWallDistanceExtrapolate(){
  Framework::DataHandle <bool> nodeisAD = socket_nodeisAD.getDataHandle();
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  Common::SafePtr<Framework::TopologicalRegionSet> cells =
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  for (CFuint iCell=0; iCell<nbCells; ++iCell){
    geoData.idx = iCell;
    CFuint nbNodesInSideRegion = 0;
    bool foundCell = false;
    const CFuint nbNodesInCell = cells->getNbNodesInGeo(iCell);
    for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
      CFuint nodeIDinC=cells->getNodeID(iCell, iNodeC);
      if ( (nodeisAD[nodeIDinC] == true) ){  
	  foundCell = true;
	}
	
	if(foundCell){
	  for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
	    CFuint nodeIDinC2=cells->getNodeID(iCell, iNodeC);

	    if(nodeisAD[nodeIDinC2]){
	      nbNodesInSideRegion= nbNodesInSideRegion+1;
	      
	    }
	  }
	  if(nbNodesInSideRegion >= 3){
	    for (CFuint iNodeC = 0; iNodeC < nbNodesInCell; ++iNodeC){
	      CFuint nodeIDinC3=cells->getNodeID(iCell, iNodeC);
		nodeisAD[nodeIDinC3] = true;
	    }
	  }
	}
      }
  }
    m_geoBuilder.releaseGE();
}


//////////////////////////////////////////////////////////////////////////////
void MeshFittingAlgorithm::updateNodePositions () {
  CFAUTOTRACE;
  CFLog(VERBOSE, "MeshFittingAlgorithm::updateNodePositions()\n"); 
  Framework::DataHandle<CFreal> nodeDistance = socket_nodeDistance.getDataHandle(); 
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle < CFreal > rhs = socket_rhs.getDataHandle();
  const CFuint totalNbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  Framework::DataHandle<CFreal> iradius = socket_iradius.getDataHandle();
  Framework::DataHandle<CFreal> skewness = socket_skewness.getDataHandle();
  Framework::DataHandle<CFreal> AR = socket_AR.getDataHandle();
  Framework::DataHandle<CFreal> isphere = socket_isphere.getDataHandle();
  const CFuint nbDims = Framework::PhysicalModelStack::getActive()->getDim();
  for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) { 
   CFreal f=1.;
    if (nodes[iNode]->isParUpdatable()) {
      Framework::Node& currNode = *nodes[iNode];
      bool exit = false;
      if (m_smoothNodalDisp){
	// FB : test case dependent

	if( nodeDistance[nodes[iNode]->getLocalID()] < 0.01 &&  insideRegion(nodes[iNode])==false ){ 				
	  f=(1./(0.001-m_acceptableDistance)) * nodeDistance[nodes[iNode]->getLocalID()] - (m_acceptableDistance/(0.001-m_acceptableDistance)); 
	}
      }
      for(CFuint iDim = 0; iDim < nbDims; ++iDim) {
	currNode[XX+iDim] =  currNode[XX+iDim]*(1.-m_meshAcceleration*f) + rhs[iNode*totalNbEqs+XX+iDim]*m_meshAcceleration*f;
      }
      
    }
  }

  /////////2d triangular 
  if(m_MQIvalue==2){
    CFuint nbPairsNodeNode = 300000;
    typedef CFMultiMap<CFuint, CFreal> MapNodeRadiusF;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeRadiusF(nbPairsNodeNode);
    typedef CFMultiMap<CFuint, CFreal> MapNodeNStateF;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNStateF(nbPairsNodeNode);
    
    Common::SafePtr<Framework::TopologicalRegionSet> cells =
      Framework::MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    geoData.trs = cells;
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
      geoData.idx = iCell;
      Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
      Framework::Node * thirdNode;
      const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
      const CFuint nbFaces = facesInCell.size(); 
      RealVector res(3);res=0.;
      std::vector<Framework::Node*>& faceNodes1 = *facesInCell[0]->getNodes();
      bool found=false;
      RealVector v1(3); v1=0.;
      RealVector v2(3); v2=0.;
      RealVector v3(3); v3=0.;
      CFreal v1state;
      if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
	v1state= (*states[iCell])[m_monitorVarID];
      }
      else{
	cf_assert(m_monitorPhysVarID < m_pdata.size());
	m_state->copyData(*states[iCell]);
	getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
	v1state  = m_pdata[m_monitorPhysVarID];
      }
      for (CFuint iFace=1; iFace<nbFaces; ++iFace){
	while(found==false){
	  std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	  if(faceNodes[0]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[0]->getLocalID() != faceNodes1[1]->getLocalID()){
	    found=true;
	    thirdNode =faceNodes[0];
	  }
	  if(faceNodes[1]->getLocalID() != faceNodes1[0]->getLocalID() && faceNodes[1]->getLocalID() != faceNodes1[1]->getLocalID()){
	    found=true;
	    thirdNode =faceNodes[1];
	  }
	}
      }
      for (CFuint iDim=0; iDim<nbDims ;++iDim){
	v1[iDim] =(*faceNodes1[0])[iDim]-(*faceNodes1[1])[iDim];
	v2[iDim] =(*faceNodes1[0])[iDim]-(*thirdNode)[iDim];
	v3[iDim]=(*faceNodes1[1])[iDim]-(*thirdNode)[iDim];
      }
      MathTools::MathFunctions::crossProd(v1, v2,res);
      
      CFreal area =res.norm2();
      CFreal perim = v1.norm2()+v2.norm2()+v3.norm2();
      CFreal radius =area/perim;
      
      
      m_mapNodeRadiusF.insert(faceNodes1[0]->getLocalID(),radius);
      m_mapNodeNStateF.insert(faceNodes1[0]->getLocalID(),v1state);
      m_mapNodeRadiusF.insert(faceNodes1[1]->getLocalID(),radius);
      m_mapNodeNStateF.insert(faceNodes1[1]->getLocalID(),v1state);
      m_mapNodeRadiusF.insert(thirdNode->getLocalID(),radius);
      m_mapNodeNStateF.insert(thirdNode->getLocalID(),v1state);
      
      
      m_geoBuilder.releaseGE();
      m_mapNodeRadiusF.sortKeys();
      m_mapNodeNStateF.sortKeys();
      
    }
    typedef MapNodeRadiusF::MapIterator mapItRf;
    typedef MapNodeNStateF::MapIterator mapItSi;
    
    for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
      if (nodes[iNode]->isParUpdatable()) {
	bool foundRi = false;
	bool foundRf = false;
	bool foundSi = false;
	bool foundSf = false;
	std::pair<mapItRf,mapItRf > itri=m_mapNodeRadius1.find(nodes[iNode]->getLocalID(), foundRi);
	std::pair<mapItRf,mapItRf > itrf=m_mapNodeRadiusF.find(nodes[iNode]->getLocalID(), foundRf);
	std::pair<mapItSi,mapItSi > itsi=m_mapNodeNState1.find(nodes[iNode]->getLocalID(), foundSi);
	std::pair<mapItSi,mapItSi > itsf=m_mapNodeNStateF.find(nodes[iNode]->getLocalID(), foundSf);
	
	cf_assert(foundRi);
	cf_assert(foundRf);
	cf_assert(foundSi);
	cf_assert(foundSf);
	
	CFreal sumRadf=0.;
	CFreal sumRadi=0.;
	CFreal sumOfOldStates=0.;
	CFreal sumOfNewStates=0.;
	
	for (mapItRf it = itrf.first; it != itrf.second; ++it){
	  sumRadf+=it->second;
	  //sizeOfSecond+=1.;
	}
	for (mapItSi it2 = itsi.first; it2 != itsi.second; ++it2){
	  sumOfOldStates+=it2->second;
	}
	for (mapItSi it3 = itsf.first; it3 != itsf.second; ++it3){
	  sumOfNewStates+=it3->second;
	  //sizeOfSecond+=1.;
	}
	for (mapItRf it1 = itri.first; it1 != itri.second; ++it1){
	  sumRadi+=it1->second;
	}
	
	CFreal finalRadius =sumRadf/sumRadi;       // not used sizeOfSecond since the connectivity does not change
	CFreal finalState = sumOfNewStates/sumOfOldStates;
       	CFreal value = (finalRadius*finalState);
	iradius[nodes[iNode]->getLocalID()] =value;	
      }
    }
  }
  ///////////////////////
   //2D quads skewness
  if(m_MQIvalue==4){
    Common::SafePtr<Framework::TopologicalRegionSet> cells =
      Framework::MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
    typedef MapNodeNode::MapIterator mapIt;
    typedef MapNodeNode::MapIterator mapItN;
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    typedef CFMultiMap<CFuint, CFreal> MapNodeSkewF;
    CFuint nbPairsNodeSkew =200000 ; 
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeSkewF(nbPairsNodeSkew);
    typedef CFMultiMap<CFuint, CFreal> MapNodeNStateF;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNStateF(nbPairsNodeSkew);
    geoData.trs = cells;
    
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
      geoData.idx = iCell;
      Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
      Framework::Node * firstNode;
      Framework::Node * secondNode;
      Framework::Node * thirdNode;
      Framework::Node * fourthNode;
      
      const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
      const CFuint nbFaces = facesInCell.size(); 
      std::vector<Framework::Node*>& faceNodes12 = *facesInCell[0]->getNodes();
      firstNode = faceNodes12[0];
      secondNode = faceNodes12[1];
      bool foundFN =false;
      bool foundSN =false;
      
      std::pair<mapIt,mapIt > itFirstNode=m_mapNodeNode1.find(firstNode->getLocalID(), foundFN);
      std::pair<mapItN,mapItN > itSecondNode=m_mapNodeNode1.find(secondNode->getLocalID(), foundSN);
      cf_assert(foundFN);
      cf_assert(foundSN);
      
      for (CFuint iFace=1; iFace<nbFaces; ++iFace){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	const CFuint nbNodesinF = faceNodes.size();
	for (mapIt itFN = itFirstNode.first; itFN != itFirstNode.second; ++itFN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itFN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	      fourthNode =  faceNodes[iNode];
	    }
	  }
	}      
	for (mapItN itSN = itSecondNode.first; itSN != itSecondNode.second; ++itSN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itSN->second == faceNodes[iNode]->getLocalID()  && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID() ){
	      thirdNode =  faceNodes[iNode];
	    }
	  }
	}
      }
      
      CFreal v1state;
      if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
	v1state= (*states[iCell])[m_monitorVarID];
      }
      else{
	cf_assert(m_monitorPhysVarID < m_pdata.size());
	m_state->copyData(*states[iCell]);
	getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
	v1state  = m_pdata[m_monitorPhysVarID];
      }
      CFreal skew;
      skew =computeSkewness2dQuads(firstNode,secondNode,thirdNode,fourthNode);
      m_mapNodeSkewF.insert(firstNode->getLocalID(),skew);
      m_mapNodeSkewF.insert(secondNode->getLocalID(),skew);
      m_mapNodeSkewF.insert(thirdNode->getLocalID(),skew);
      m_mapNodeSkewF.insert(fourthNode->getLocalID(),skew);
      m_mapNodeNStateF.insert(firstNode->getLocalID(),v1state);
      m_mapNodeNStateF.insert(secondNode->getLocalID(),v1state);
      m_mapNodeNStateF.insert(thirdNode->getLocalID(),v1state);
      m_mapNodeNStateF.insert(fourthNode->getLocalID(),v1state);
      m_mapNodeSkewF.sortKeys();
      m_mapNodeNStateF.sortKeys();
      
      m_geoBuilder.releaseGE();
      
    }

    
    typedef MapNodeSkewF::MapIterator mapItSkf;
    typedef MapNodeNStateF::MapIterator mapItSi;
    
    for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
     if (nodes[iNode]->isParUpdatable() ){  //&& isNodeLocked(nodes[iNode]) == false){
       bool foundRi = false;
       bool foundRf = false;
       bool foundSi = false;
       bool foundSf = false;
       
       std::pair<mapItSkf,mapItSkf > itrf=m_mapNodeSkewF.find(nodes[iNode]->getLocalID(), foundRf);
       std::pair<mapItSi,mapItSi > itsi=m_mapNodeNState1.find(nodes[iNode]->getLocalID(), foundSi);
       std::pair<mapItSi,mapItSi > itsf=m_mapNodeNStateF.find(nodes[iNode]->getLocalID(), foundSf);
       std::pair<mapItSkf,mapItSkf > itri=m_mapNodeSkew1.find(nodes[iNode]->getLocalID(), foundRi);
       cf_assert(foundRf);
       cf_assert(foundSi);
       cf_assert(foundSf);
       cf_assert(foundRi);
       CFreal sumAngf=0.;
       CFreal sumAngi=0.;
       CFreal sumOfOldStates=0.;
       CFreal sumOfNewStates=0.;
       CFreal sizeOfSecond=0.;
       CFreal sizeOfSecond1=0.;
       
       for (mapItSkf it = itrf.first; it != itrf.second; ++it){
	 sumAngf+=it->second;
	 
	 sizeOfSecond +=1.;
       }
       for (mapItSi it2 = itsi.first; it2 != itsi.second; ++it2){
	 sumOfOldStates+=it2->second;
       }
       for (mapItSi it3 = itsf.first; it3 != itsf.second; ++it3){
	 sumOfNewStates+=it3->second;
	 sizeOfSecond1 +=1.;
       }
       for (mapItSkf it1 = itri.first; it1 != itri.second; ++it1){
	 sumAngi+=it1->second;
       }
       
       CFreal finalSkew =std::abs(sumAngf/sizeOfSecond-sumAngi/sizeOfSecond1);              
       CFreal finalState = sumOfNewStates/sumOfOldStates;
       
       CFreal value = (finalSkew*finalState);
       skewness[nodes[iNode]->getLocalID()] =value;	
     }
    }
    
  }


 ///////////////////////
   //2D quads AR
  if(m_MQIvalue==3){
    Common::SafePtr<Framework::TopologicalRegionSet> cells =
      Framework::MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    typedef CFMultiMap<CFuint, CFuint> MapNodeNode;
    typedef MapNodeNode::MapIterator mapIt;
    typedef MapNodeNode::MapIterator mapItN;
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    typedef CFMultiMap<CFuint, CFreal> MapNodeARF;
    CFuint nbPairsNodeSkew =200000 ; 
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeARF(nbPairsNodeSkew);
    typedef CFMultiMap<CFuint, CFreal> MapNodeNStateF;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNStateF(nbPairsNodeSkew);
    geoData.trs = cells;
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
      geoData.idx = iCell;
      Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
      Framework::Node * firstNode;
      Framework::Node * secondNode;
      Framework::Node * thirdNode;
      Framework::Node * fourthNode;
      
      const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
      const CFuint nbFaces = facesInCell.size(); 
      std::vector<Framework::Node*>& faceNodes12 = *facesInCell[0]->getNodes();
      firstNode = faceNodes12[0];
      secondNode = faceNodes12[1];
      bool foundFN =false;
      bool foundSN =false;
      
      std::pair<mapIt,mapIt > itFirstNode=m_mapNodeNode1.find(firstNode->getLocalID(), foundFN);
      std::pair<mapItN,mapItN > itSecondNode=m_mapNodeNode1.find(secondNode->getLocalID(), foundSN);
      cf_assert(foundFN);
      cf_assert(foundSN);
      
      for (CFuint iFace=1; iFace<nbFaces; ++iFace){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	const CFuint nbNodesinF = faceNodes.size();
	for (mapIt itFN = itFirstNode.first; itFN != itFirstNode.second; ++itFN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itFN->second == faceNodes[iNode]->getLocalID() && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID()){
	      fourthNode =  faceNodes[iNode];
	    }
	  }
	}      
	for (mapItN itSN = itSecondNode.first; itSN != itSecondNode.second; ++itSN) {
	  for(CFuint iNode=0; iNode<nbNodesinF ; ++iNode){
	    if(itSN->second == faceNodes[iNode]->getLocalID()  && faceNodes[iNode]->getLocalID()!=firstNode->getLocalID()  && faceNodes[iNode]->getLocalID()!=secondNode->getLocalID() ){
	      thirdNode =  faceNodes[iNode];
	    }
	  }
	}
      }
      
      CFreal v1state;
      if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
	v1state= (*states[iCell])[m_monitorVarID];
      }
      else{
	cf_assert(m_monitorPhysVarID < m_pdata.size());
	m_state->copyData(*states[iCell]);
	getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
	v1state  = m_pdata[m_monitorPhysVarID];
      }
      CFreal AR;
      AR =computeAspectRatio2dQuads(firstNode,secondNode,thirdNode,fourthNode);
      m_mapNodeARF.insert(firstNode->getLocalID(),AR);
      m_mapNodeARF.insert(secondNode->getLocalID(),AR);
      m_mapNodeARF.insert(thirdNode->getLocalID(),AR);
      m_mapNodeARF.insert(fourthNode->getLocalID(),AR);
      m_mapNodeNStateF.insert(firstNode->getLocalID(),v1state);
      m_mapNodeNStateF.insert(secondNode->getLocalID(),v1state);
      m_mapNodeNStateF.insert(thirdNode->getLocalID(),v1state);
      m_mapNodeNStateF.insert(fourthNode->getLocalID(),v1state);
      m_mapNodeARF.sortKeys();
      m_mapNodeNStateF.sortKeys();
      
      m_geoBuilder.releaseGE();
      
    }
    
    
    typedef MapNodeARF::MapIterator mapItSkf;
    typedef MapNodeNStateF::MapIterator mapItSi;
    
    for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
      if (nodes[iNode]->isParUpdatable() ){ 
	bool foundRi = false;
	bool foundRf = false;
	bool foundSi = false;
	bool foundSf = false;
	
	std::pair<mapItSkf,mapItSkf > itrf=m_mapNodeARF.find(nodes[iNode]->getLocalID(), foundRf);
	std::pair<mapItSi,mapItSi > itsi=m_mapNodeNState1.find(nodes[iNode]->getLocalID(), foundSi);
	std::pair<mapItSi,mapItSi > itsf=m_mapNodeNStateF.find(nodes[iNode]->getLocalID(), foundSf);
	std::pair<mapItSkf,mapItSkf > itri=m_mapNodeAR1.find(nodes[iNode]->getLocalID(), foundRi);
	cf_assert(foundRf);
	cf_assert(foundSi);
	cf_assert(foundSf);
	cf_assert(foundRi);
	CFreal sumARf=0.;
	CFreal sumARi=0.;
	CFreal sumOfOldStates=0.;
	CFreal sumOfNewStates=0.;
	//CFreal sizeOfSecond=0.;
	//CFreal sizeOfSecond1=0.;
	
	for (mapItSkf it = itrf.first; it != itrf.second; ++it){
	  sumARf+=it->second;
	  //sizeOfSecond +=1.;
	}
	for (mapItSi it2 = itsi.first; it2 != itsi.second; ++it2){
	  sumOfOldStates+=it2->second;
	}
	for (mapItSi it3 = itsf.first; it3 != itsf.second; ++it3){
	  sumOfNewStates+=it3->second;
	  //sizeOfSecond1 +=1.;
	}
	for (mapItSkf it1 = itri.first; it1 != itri.second; ++it1){
	  sumARi+=it1->second;
	}
	
	CFreal finalAR =(sumARi/sumARf);       // not used sizeOfSecond since the connectivity does not change
	CFreal finalState = sumOfNewStates/sumOfOldStates;
	CFreal value = (finalAR*finalState);
	AR[nodes[iNode]->getLocalID()] =value;	
      }
    }
  }
//3D tetrahedral
  if(m_MQIvalue == 5 ){
    Common::SafePtr<Framework::TopologicalRegionSet> cells = 
      Framework::MeshDataStack::getActive()->getTrs("InnerCells");
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    Framework::CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    geoData.trs = cells;
    CFuint nbPairsNodeNode1 =300000; 
    
    typedef CFMultiMap<CFuint, CFuint> MapCellNode;
    typedef CFMultiMap<CFuint, CFreal > MapNodeNState;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeNStateF(nbPairsNodeNode1);
    typedef CFMultiMap<CFuint, CFreal > MapNodeRS;
    Common::CFMultiMap<CFuint,CFreal>  m_mapNodeTSF(nbPairsNodeNode1);
    typedef MapCellNode::MapIterator mapItc;
    for (CFuint iCell=0; iCell<nbCells; ++iCell){
      geoData.idx = iCell;
      std::vector<CFuint> nodeIDs;
      nodeIDs.clear();
      CFreal facesArea = 0.;
      CFreal volume =0.;
      CFreal radius = 0.;
      Framework::GeometricEntity *const currCell = m_geoBuilder.buildGE();
      bool foundC = false;
      std::pair<mapItc,mapItc > iteC=m_mapCellNode1.find(iCell, foundC);
      cf_assert(foundC);
      for (mapItc it2 = iteC.first; it2 != iteC.second; ++it2){
	nodeIDs.push_back(it2->second);
      }
      volume=ComputeTvolume(nodeIDs[0],nodeIDs[1],nodeIDs[2],nodeIDs[3]);
      const std::vector<Framework::GeometricEntity*>& facesInCell = *currCell->getNeighborGeos();
      const CFuint nbFaces = facesInCell.size(); 
      for (CFuint iFace=1; iFace<nbFaces; ++iFace){
	std::vector<Framework::Node*>& faceNodes = *facesInCell[iFace]->getNodes();
	Framework::Node* node1 = faceNodes[0];
	Framework::Node* node2 = faceNodes[1];
	Framework::Node* node3 = faceNodes[2];
	facesArea+= ComputeTFacesurface(node1,node2,node3);
      }
      radius = 3.*volume/facesArea;
      CFreal v1state ;
      if (m_monitorPhysVarID == std::numeric_limits<CFuint>::max()) {
	v1state= (*states[iCell])[m_monitorVarID];
      }
      else{
	cf_assert(m_monitorPhysVarID < m_pdata.size());
	m_state->copyData(*states[iCell]);
	getMethodData().getUpdateVarSet()->computePhysicalData(*m_state, m_pdata);
	v1state  = m_pdata[m_monitorPhysVarID];
    }
      for(CFuint i=0; i<nodeIDs.size() ; ++i){
	m_mapNodeTSF.insert(nodeIDs[i],radius);
	m_mapNodeNStateF.insert(nodeIDs[i],v1state);
      }
      m_mapNodeTSF.sortKeys();
      m_mapNodeNStateF.sortKeys();
      
      m_geoBuilder.releaseGE();
      
    }
    typedef MapNodeRS::MapIterator mapItRadf;
    typedef MapNodeNState::MapIterator mapItSi;
    for (CFuint iNode = 0; iNode < nodes.size(); ++iNode) {  
      if (nodes[iNode]->isParUpdatable() ){
	bool foundRi = false;
	bool foundRf = false;
	bool foundSi = false;
	bool foundSf = false;
	std::pair<mapItRadf,mapItRadf > itrf=m_mapNodeTSF.find(nodes[iNode]->getLocalID(), foundRf);
	std::pair<mapItSi,mapItSi > itsi=m_mapNodeNState1.find(nodes[iNode]->getLocalID(), foundSi);
	std::pair<mapItSi,mapItSi > itsf=m_mapNodeNStateF.find(nodes[iNode]->getLocalID(), foundSf);
	std::pair<mapItRadf,mapItRadf > itri=m_mapNodeTS1.find(nodes[iNode]->getLocalID(), foundRi);
	cf_assert(foundRf);
	cf_assert(foundSi);
	cf_assert(foundSf);
	cf_assert(foundRi);
	CFreal sumRadf=0.;
	CFreal sumRadi=0.;
	CFreal sumOfOldStates=0.;
	CFreal sumOfNewStates=0.;
	//CFreal sizeOfSecond=0.;
	//CFreal sizeOfSecond1=0.;
	for (mapItRadf it = itrf.first; it != itrf.second; ++it){
	  sumRadf+=it->second;
	  //sizeOfSecond +=1.;
	}
	for (mapItSi it2 = itsi.first; it2 != itsi.second; ++it2){
	  sumOfOldStates+=it2->second;
	}
	for (mapItSi it3 = itsf.first; it3 != itsf.second; ++it3){
	  sumOfNewStates+=it3->second;
	  //sizeOfSecond1 +=1.;
	}
	for (mapItRadf it1 = itri.first; it1 != itri.second; ++it1){
	  sumRadi+=it1->second;
	}
	
	CFreal finalRadius =(sumRadf/sumRadi);      
	CFreal finalState = sumOfNewStates/sumOfOldStates;
	
	CFreal value = (finalRadius*finalState);
	isphere[nodes[iNode]->getLocalID()] =value;	
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
