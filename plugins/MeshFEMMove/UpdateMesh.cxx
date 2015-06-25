#include "Environment/FileHandlerOutput.hh"
#include "Common/ProcessInfo.hh"
#include "Common/OSystem.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/Node.hh"
#include "Framework/PathAppender.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/State.hh"

#include "MeshFEMMove/MeshFEMMove.hh"
#include "MeshFEMMove/UpdateMesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateMesh, FEMMoveData, MeshFEMMoveModule> updateMeshProvider("UpdateMesh");

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption<CFint>("ScreenOutputPrecision","Screen Output Precision.");
}

//////////////////////////////////////////////////////////////////////////////

UpdateMesh::UpdateMesh(const std::string& name) :
    FEMMoveCom(name),
    socket_nodes("nodes"),
    socket_states("states"),
    socket_otherNodes("nodes")
{

   addConfigOptionsTo(this);

  _precision = 8;
  setParameter("ScreenOutputPrecision",&_precision);

}

//////////////////////////////////////////////////////////////////////////////

UpdateMesh::~UpdateMesh()
{
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FEMMoveCom::configure(args);

  socket_otherNodes.setDataSocketNamespace(getMethodData().getOtherNamespace());
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateMesh::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_otherNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::execute()
{
  CFAUTOTRACE;

  MultiMethodHandle<ConvergenceMethod> convergenceMethod = getMethodData().getConvergenceMethod();
  cf_assert(convergenceMethod.size() == 1);

  _stopwatch.start();

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  // reset to 0 the residual
  subSysStatus->resetResidual();

  // set the handle to the global states and nodes in the
  // convergence method
///@todo here or in standard subsystem...change this!!!
  //_convergenceMethod[0]->setGlobalStates(statedata);***************************
  //_convergenceMethod[0]->setGlobalNodes(nodedata);***************************

  convergenceMethod[0]->takeStep();

  writeOnScreen();

  moveMeshNodes();
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::moveMeshNodes()
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<Node*, GLOBAL> otherNodes = socket_otherNodes.getDataHandle();

  cf_assert(nodes.size() == states.size());
  cf_assert(nodes.size() == otherNodes.size());

  std::string name = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> otherPhysModel = PhysicalModelStack::getInstance().getEntryByNamespace(otherNsp);
  
  cf_assert(otherPhysModel->getDim() == PhysicalModelStack::getActive()->getNbEq());

  for(CFuint iNode = 0; iNode < otherNodes.size(); ++iNode)
  {
     *(otherNodes[iNode]) += *(states[iNode]);
  }

  if(SubSystemStatusStack::getActive()->isSubIterationLastStep())
  {
    for(CFuint iNode = 0; iNode < nodes.size(); ++iNode)
    {
      *(nodes[iNode]) += *(states[iNode]);
    }
  }

    for(CFuint iState = 0; iState < states.size(); ++iState)
    {
      *(states[iState]) = 0.;
    }

//   for(CFuint iState = 0; iState < states.size(); ++iState)
//   {
//     if (states[iState]->isParUpdatable())
//       *(states[iState]) = 1.;
//     if (!states[iState]->isParUpdatable())
//       *(states[iState]) = 0.;
//   }

}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::writeOnScreen()
{
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  ostringstream out;

  // Print number of iterations
  out << "FEM Mesh Movement - Iter: ";
  out.width(5);
  out <<  subSysStatus->getNbIter() << " ";

  // Print Residual
  out << " - Res:" << " ";
  out.width(3 + _precision);
  out.precision(_precision);
  out << subSysStatus->getResidual();

  // Print CPU Time
  out << " -  CPUTime: ";
  out.width(3);
  out << _stopwatch.read() << " ";
  out.width(3);

  // Print Memory Usage
  out << " - Mem: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n";

  CFout << out.str() << CFendl;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD
