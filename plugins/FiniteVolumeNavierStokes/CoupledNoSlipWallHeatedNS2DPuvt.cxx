#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/CoupledNoSlipWallHeatedNS2DPuvt.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledNoSlipWallHeatedNS2DPuvt, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
coupledNoSlipWallHeatedNS2DPuvtFVMCCProvider("CoupledNoSlipWallHeatedNS2DPuvtFVMCC");

//////////////////////////////////////////////////////////////////////////////

void CoupledNoSlipWallHeatedNS2DPuvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Interface","Name of the Interface.");
  options.addConfigOption< CFuint >("DefaultIterations","Number of iterations during which to use the default values.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledNoSlipWallHeatedNS2DPuvt::CoupledNoSlipWallHeatedNS2DPuvt(const std::string& name) :
  NoSlipWallHeatedNSPvt(name),
  _sockets()
{
   addConfigOptionsTo(this);

  _interfaceName = "";
  setParameter("Interface",&_interfaceName);

  _defaultIterations = 0;
  setParameter("DefaultIterations",&_defaultIterations);

  _setIndex = false;
}

//////////////////////////////////////////////////////////////////////////////

CoupledNoSlipWallHeatedNS2DPuvt::~CoupledNoSlipWallHeatedNS2DPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledNoSlipWallHeatedNS2DPuvt::configure ( Config::ConfigArgs& args )
{
  NoSlipWallHeatedNSPvt::configure(args);

  const std::string nameSpace = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string> trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];

    const std::string baseSocketName =
      "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";
    
    _sockets.createSocketSink<CFreal>(baseSocketName + "ISACCEPTED");
    _sockets.createSocketSink<RealVector>(baseSocketName + "DATA");
    _sockets.createSocketSink<RealVector>(baseSocketName + "COORD");
  }
}

//////////////////////////////////////////////////////////////////////////////

void CoupledNoSlipWallHeatedNS2DPuvt::setup()
{
  NoSlipWallHeatedNSPvt::setup();
  
  const CFuint nbDim = Framework::PhysicalModelStack::getActive()->getDim();
  m_midNode.resize(nbDim);

  //Create a Map to get back the index of the node in the TRS list from its LocalID
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();

  const CFuint nbNodesInTRS = trsNodes->size();
  for (CFuint iNode = 0; iNode < nbNodesInTRS; ++iNode) {
    const CFuint nodeID = (*trsNodes)[iNode];
    _trsNodeIDMap.insert(nodeID, iNode);
  }
  _trsNodeIDMap.sortKeys();



}

//////////////////////////////////////////////////////////////////////////////

void CoupledNoSlipWallHeatedNS2DPuvt::setIndex()
{
  CFAUTOTRACE;

  ///@todo move this to setup
  Common::SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  const std::string trsName = trs->getName();
  const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
  const std::string nameSpace = getMethodData().getNamespace();

  const std::string baseSocketName =
    "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";

  DataHandle<CFreal> isAccepted =
    _sockets.getSocketSink<CFreal>(baseSocketName + "ISACCEPTED")->getDataHandle();

  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();
  const CFuint nbNodesInTRS = trsNodes->size();

  cf_assert(isAccepted.size() == trsNodes->size());

  _coupledDataID.resize(nbNodesInTRS);
  CFuint idx = 0;
  for(CFuint iNode=0; iNode < nbNodesInTRS; ++iNode)
  {
    if(isAccepted[iNode] >= 0.)
    {
      _coupledDataID[iNode] = idx;
      idx++;
    }
    else{
      _coupledDataID[iNode] = -1;
    }
  }

  _setIndex = true;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
CoupledNoSlipWallHeatedNS2DPuvt::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = NoSlipWallHeatedNSPvt::needsSockets();
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result2 = _sockets.getAllSinkSockets();

  for(CFuint i=0;i<result2.size();++i)
  {
    result.push_back(result2[i]);
  }

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void CoupledNoSlipWallHeatedNS2DPuvt::setGhostState(GeometricEntity *const face)
{

  if(!_setIndex) setIndex();

  Common::SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();

  const std::string trsName = trs->getName();
  const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
  const std::string nameSpace = getMethodData().getNamespace();

  const std::string baseSocketName =
    "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";

  DataHandle<RealVector> interfaceData =
    _sockets.getSocketSink<RealVector>(baseSocketName + "DATA")->getDataHandle();
  DataHandle<CFreal> isAccepted =
    _sockets.getSocketSink<CFreal>(baseSocketName + "ISACCEPTED")->getDataHandle();
  DataHandle<RealVector> coordinates =
    _sockets.getSocketSink<RealVector>(baseSocketName + "COORD")->getDataHandle();

  ///Check that the datahandle has same size as the TRS Nodes
  cf_assert(isAccepted.size() == trsNodes->size());

  vector<Node*> *const faceNodes = face->getNodes();
  const CFuint nbFaceNodes = faceNodes->size();

  //Get the transfered values at the nodes of the face
  //First find the corresponding index
  vector<CFuint> faceNodesTrsID(nbFaceNodes);
  for(CFuint iNode=0; iNode < nbFaceNodes; ++iNode)
  {
    faceNodesTrsID[iNode] = _trsNodeIDMap.find((*faceNodes)[iNode]->getLocalID());
  }
  // then get the values at the face nodes from the coupled DataHandle
  RealVector faceNodesFlux(nbFaceNodes);

  for(CFuint iNode=0; iNode < nbFaceNodes; ++iNode)
  {
    if((isAccepted[faceNodesTrsID[iNode]] >= 0.) &&
       (SubSystemStatusStack::getActive()->getNbIter() > _defaultIterations))
    {
      cf_assert(interfaceData[_coupledDataID[faceNodesTrsID[iNode]]].size() == 2);
      cf_assert(interfaceData[_coupledDataID[faceNodesTrsID[iNode]]][0] > 0.);

      faceNodesFlux[iNode] = -1.*interfaceData[_coupledDataID[faceNodesTrsID[iNode]]][1];
    }
    else
    {
      //set the default value
      faceNodesFlux[iNode] = _heatFlux;
    }
  }


  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // this middle node is by construction on the boundary face
  m_midNode = 0.5*(innerState->getCoordinates() + ghostState->getCoordinates());

  //Compute the wall value (at m_midNode)
  CFreal wallFlux = 0.;
  CFreal totalDistance = 0.;
  for(CFuint iNode=0; iNode < nbFaceNodes; ++iNode)
  {
    const CFreal distance =
      MathFunctions::getDistance(m_midNode, *(face->getNode(iNode)));
    wallFlux += (faceNodesFlux[iNode] * distance);
    totalDistance += distance;
  }
  wallFlux /= totalDistance;


  ///Set Tghost such that Q = lambda * (Tghost - Tin)/(distance_G_In)
  ///-> set TG = Q*d/lambda + Tin
  const CFreal distance = MathFunctions::getDistance(innerState->getCoordinates(),ghostState->getCoordinates());
  _diffusiveVarSet->setWallDistance(distance/2.);
  const CFreal dynamicViscosity = _diffusiveVarSet->getDynViscosity(*innerState, _dummyGradients);

  const CFreal lambda = _diffusiveVarSet->getThermConductivity(*innerState,dynamicViscosity);
// unused //  const CFreal area = MathFunctions::getDistance(*(face->getNode(0)), *(face->getNode(1)));
  (*ghostState)[0] = (*innerState)[0];
  (*ghostState)[1] = 2.*_xWallVelocity - (*innerState)[1];
  (*ghostState)[2] = 2.*_yWallVelocity - (*innerState)[2];
  (*ghostState)[3] = (wallFlux * distance/ lambda) + (*innerState)[3];

/*  CFout << "distance: " << distance << "\n";
  CFout << "wallFlux: " << wallFlux << "\n";
  CFout << "area: " << area << "\n";
  CFout << "lambda: " << lambda << "\n";
  CFout << "Inner T: " << (*innerState)[3] << "\n";
  CFout << "Ghost T: " << (*ghostState)[3] << "\n";*/
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
