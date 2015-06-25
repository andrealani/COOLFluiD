#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"
#include "CoupledPrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/SelfRegistPtr.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledPrepare, SpringAnalogyData, MeshAdapterSpringAnalogyModule> CoupledPrepareProvider("CoupledPrepare");

//////////////////////////////////////////////////////////////////////////////

void CoupledPrepare::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::string >("DataType","Type of Data Received Coordinates or Displacements.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledPrepare::CoupledPrepare(const std::string& name) :
SpringAnalogyCom(name),
  _sockets(),
  socket_nodes("nodes"),
  socket_nodalDisplacements("nodalDisplacements")
{
   addConfigOptionsTo(this);

  _interfaceName = "";
   setParameter("Interface",&_interfaceName);

  _dataTypeReceived = "";
   setParameter("DataType",&_dataTypeReceived);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledPrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_nodalDisplacements);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CoupledPrepare::configure ( Config::ConfigArgs& args )
{
  SpringAnalogyCom::configure(args);

  const std::string nameSpace = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
  
  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];
    std::string socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_ISACCEPTED";
    _sockets.createSocketSink<CFreal>(socketName);

    socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_DATA";
    _sockets.createSocketSink<RealVector>(socketName);
  }

}

//////////////////////////////////////////////////////////////////////////////

void CoupledPrepare::execute()
{
  CFAUTOTRACE;

  moveBoundaries();

//  computeWorstQualityCells();

 }

//////////////////////////////////////////////////////////////////////////////

void CoupledPrepare::moveBoundaries()
{
  CFAUTOTRACE;

  /// Do the coupling part: move the boundary nodes by getting the datahandle
  const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> displacements = socket_nodalDisplacements.getDataHandle();

  vector< Common::SafePtr<TopologicalRegionSet> > trs = getTrsList();
  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    const std::string currentTrsName = getTrsName(iTRS);

    ///Read the file and put it into the datahandle
    const std::string trsName = getTrsName(iTRS);
    const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
    const std::string nameSpace = getMethodData().getNamespace();

    std::string socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_DATA";
    DataHandle<RealVector> interfaceData = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();

    socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_ISACCEPTED";
    DataHandle<CFreal> isAccepted = _sockets.getSocketSink<CFreal>(socketName)->getDataHandle();

    // Modify the coordinates with the new coordinates received
    SafePtr<TopologicalRegionSet> currentTrs = getCurrentTRS();

    SafePtr<vector<CFuint> > trsNodes = currentTrs->getNodesInTrs();
    vector<CFuint>::iterator itd;

    CFuint currentStep = getMethodData().getCurrentStep();
    CFuint totalSteps = getMethodData().getTotalNbSteps();

    RealVector diffCoord(PhysicalModelStack::getActive()->getDim());

    CFuint idx = 0;
    CFuint totalIdx = 0;
    for (itd = trsNodes->begin(); itd != trsNodes->end(); ++itd)
      {
      const CFuint localID = *itd;
      Node* currNode = nodes[localID];

      if(isAccepted[totalIdx] >= 0.)
      {
        // Set the coordinates to the new value
        ///@todo remove, this is only to test...
        if(_dataTypeReceived == "Test")
        {
          CFreal movementCoef = 1./(totalSteps-(currentStep-1));
          diffCoord = 0.;
          diffCoord[0] = movementCoef *((0.025*interfaceData[idx][0]) - (*currNode)[0]);
          (*currNode) += diffCoord;
          displacements[localID] = diffCoord;
        }
        if(_dataTypeReceived == "Coordinates")
        {
          // Check that dimension of the Data Transfered == Dimension
          cf_assert(interfaceData[idx].size() == PhysicalModelStack::getActive()->getDim());

          CFreal movementCoef = 1/(totalSteps-(currentStep-1));
          diffCoord = movementCoef * (interfaceData[idx] - (*currNode));
          (*currNode) += diffCoord;
          displacements[localID] = diffCoord;
        }
        if(_dataTypeReceived == "Displacements")
        {
          // Check that dimension of the Data Transfered == Dimension
          cf_assert(interfaceData[idx].size() == PhysicalModelStack::getActive()->getDim());

          CFreal movementCoef = 1/totalSteps;
          diffCoord = movementCoef * interfaceData[idx];
          (*currNode) += diffCoord;
          displacements[localID] = diffCoord;
        }
        idx++ ;
      }
      totalIdx++;
    }
  }
}

// //////////////////////////////////////////////////////////////////////
//
// void CoupledPrepare::computeWorstQualityCells()
// {
//   Common::SelfRegistPtr<Framework::QualityCalculator> qualityComputer =
//     Environment::Factory<Framework::QualityCalculator>::getInstance().
//     getProvider("Concrete")->create("Concrete");
//   cf_assert(qualityComputer.isNotNull());
//
//   /// Get the datahandle of pointers to GeomEntity
//   DataHandle< std::vector<GeometricEntity*> > nodeToCellConn = socket_nodeToCellConn.getDataHandle();
//
//   DataHandle< CFreal> qualityNode = socket_qualityNode.getDataHandle();
//
//   //Loop over Nodes
//   // set the list of faces
//   std::vector< Common::SafePtr<TopologicalRegionSet> > alltrs = MeshDataStack::getActive()->getTrsList();
//   std::vector< Common::SafePtr<TopologicalRegionSet> >::iterator itrs;
//   for (itrs = alltrs.begin(); itrs != alltrs.end(); ++itrs) {
//     if ((*itrs)->getName() == "InnerCells") {
//       SafePtr<vector<CFuint> > nodesInTRS = (*itrs)->getNodesInTrs();
//       vector<CFuint>::iterator it;
//       for (it = nodesInTRS->begin(); it != nodesInTRS->end(); ++it) {
//         const CFuint localID = (*it);
//         vector<GeometricEntity*>& cells = nodeToCellConn[localID];
//         CFreal worseQuality = 0.;
//         for (CFuint iCell=0; iCell < cells.size(); ++iCell) {
//           CFreal quality = qualityComputer->computeQuality(cells[iCell]);
//           if (quality > worseQuality) worseQuality = quality;
//           if (quality < 0) worseQuality = MathTools::MathConsts::CFrealMax();
//         }
//         //cout << "--------- Worse quality : " << worseQuality << endl;
//         if(qualityNode[localID]!=0.) qualityNode[localID] = worseQuality;
//       }
//     }
//   }
//
// }

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD
