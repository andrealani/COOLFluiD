#include "FiniteVolume/FiniteVolume.hh"
#include "CoupledSuperInlet_NodalFVMCC.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledSuperInlet_NodalFVMCC,
                                 CellCenterFVMData,
                                FiniteVolumeModule>
CoupledSuperInlet_NodalFVMCCProvider("CoupledSuperInlet_NodalFVMCC");

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet_NodalFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledSuperInlet_NodalFVMCC::CoupledSuperInlet_NodalFVMCC(const std::string& name) :
  FVMCC_BC(name),
  _sockets()
{
  addConfigOptionsTo(this);

  _interfaceName = "";
  setParameter("Interface",&_interfaceName);

  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

CoupledSuperInlet_NodalFVMCC::~CoupledSuperInlet_NodalFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledSuperInlet_NodalFVMCC::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();

//  copy(_sockets.getAllSinkSockets().begin(), _sockets.getAllSinkSockets().end(), back_inserter(result));

  std::vector<Common::SafePtr<BaseDataSocketSink> > result2 = _sockets.getAllSinkSockets();

  for(CFuint i=0;i<result2.size();++i)
  {
    result.push_back(result2[i]);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet_NodalFVMCC::setup()
{
  CFAUTOTRACE;

  FVMCC_BC::setup();

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

void CoupledSuperInlet_NodalFVMCC::setGhostState(GeometricEntity *const face)
 {
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();

  const std::string trsName = trs->getName();
  const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
  const std::string nameSpace = getMethodData().getNamespace();

  const std::string baseSocketName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";

  std::string socketName = baseSocketName + "DATA";
  DataHandle<RealVector> interface = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();
  socketName = baseSocketName + "ISACCEPTED";
  DataHandle<CFreal> isAccepted = _sockets.getSocketSink<CFreal>(socketName)->getDataHandle();
  socketName = baseSocketName + "COORD";
  DataHandle<RealVector> coordinates = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();

  ///Check that the datahandle has same size as the TRS
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
  //Then get the values at the face nodes from the DataHandle
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  vector< RealVector > faceNodesValues(nbFaceNodes);
  for(CFuint iNode=0; iNode < nbFaceNodes; ++iNode)
  {
    if((isAccepted[faceNodesTrsID[iNode]] >= 0.) && (interface[faceNodesTrsID[iNode]].size() > 0))
    {
      faceNodesValues[iNode].resize(interface[faceNodesTrsID[iNode]].size());
      faceNodesValues[iNode] = interface[faceNodesTrsID[iNode]];
    }
    else
    {
      faceNodesValues[iNode].resize(nbEqs);
      _vFunction.evaluate((*((*faceNodes)[iNode])),faceNodesValues[iNode]);
    }
  }

  //Then compute the value on the face from the nodal values (average of the nodal values)
//   const CFuint transferedVectorSize = faceNodesValues[0].size();
//
//   if(transferedVectorSize == nbEqs)
//   {
//     CFout << "CoupledSuperInlet_NodalFVMCC: Transfered Vector Size different from the number of Equations" << "\n" ;
//     throw BadValueException (FromHere(),"CoupledSuperInlet_NodalFVMCC: Transfered Vector Size different from the number of Equations");
//   }

  State value;
  RealVector adimValue(nbEqs);

  value = 0.;
  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    for(CFuint iNode=0; iNode < nbFaceNodes; ++iNode)
    {
      value[iEq] += (faceNodesValues[iNode])[iEq];
    }
    value[iEq] /= nbFaceNodes;
  }

  //Adimensionalize
  getMethodData().getUpdateVar()->setAdimensionalValues(value, adimValue);

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the ghost state
  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    (*ghostState)[iEq] = 2.0*value[iEq] - (*innerState)[iEq];
  }

}

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet_NodalFVMCC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FVMCC_BC::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);

  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  const std::string nameSpace = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
  
  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];
    const std::string baseSocketName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_";

    std::string socketName = baseSocketName + "DATA";
    _sockets.createSocketSink<RealVector>(socketName);
    socketName = baseSocketName + "ISACCEPTED";
    _sockets.createSocketSink<bool>(socketName);
    socketName = baseSocketName + "COORD";
    _sockets.createSocketSink<RealVector>(socketName);
  }

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
