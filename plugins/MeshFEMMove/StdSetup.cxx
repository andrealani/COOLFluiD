#include "MeshFEMMove/MeshFEMMove.hh"
#include "StdSetup.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, FEMMoveData, MeshFEMMoveModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  FEMMoveCom(name),
  socket_femStates("states"),
  socket_femNodes("nodes"),
  socket_otherStates("states"),
  socket_otherNodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::configure ( Config::ConfigArgs& args )
{
  FEMMoveCom::configure(args);

  socket_otherStates.setDataSocketNamespace(getMethodData().getOtherNamespace());
  socket_otherNodes.setDataSocketNamespace(getMethodData().getOtherNamespace());

//   socket_femNodes.linkToSink(socket_nodes);
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  ///@todo this has to be done for the subsystem

  std::string name = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<SubSystemStatus> otherSubSystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  otherSubSystemStatus->setMovingMesh(true);

//   DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
//
//   //Create a new set of states
//   DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
//   const CFuint nbNodes = nodes.size();
//
//   states.resize(nbNodes);
//   for (CFuint i = 0; i < nbNodes; ++i) {
//     states[i] = new State();
//
// //set the nodeID into the stateID
// //set the global nodeID into the global stateID
//   }
//
//   SafePtr<TopologicalRegionSet> elements =
//     MeshDataStack::getActive()->getTrs("InnerCells");
//
//   // determining the space discretization method
//   DataHandle<State*> otherStates = socket_otherStates.getDataHandle();
//   const CFuint nbOtherStates = otherStates.size();
//   std::string spaceMethodName = (nbOtherStates == elements->getLocalNbGeoEnts()) ?
//     "FVMCC" : "FEM";
//
//   if(spaceMethodName != "FEM")
//   {
//     createExtraConnectivity();
//   }

}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::createExtraConnectivity()
{

//Copy the cellNodes connectivity table into the cellStates connectivity

// Common::SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes");
//
// ConnTable* cellStates = cellNodes;
// storeConnectivity(cellStates);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSetup::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_femNodes);
  result.push_back(&socket_femStates);
  result.push_back(&socket_otherStates);
  result.push_back(&socket_otherNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD
