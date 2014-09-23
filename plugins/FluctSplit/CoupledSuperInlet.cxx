#include "FluctSplit/FluctSplit.hh"
#include "CoupledSuperInlet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledSuperInlet, FluctuationSplitData, FluctSplitModule> coupledSuperInletProvider("CoupledSuperInlet");

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledSuperInlet::CoupledSuperInlet(const std::string& name) :
  FluctuationSplitCom(name),
  _sockets(),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_states("states")
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

CoupledSuperInlet::~CoupledSuperInlet()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledSuperInlet::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "CoupledSuperInlet::execute() called for TRS: "
  << trs->getName() << "\n");

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  Common::SafePtr< vector<CFuint> > const trsStates = trs->getStatesInTrs();

  //Build the name of the Datahandle to get data from
  const std::string trsName = trs->getName();
  const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
  const std::string nameSpace = getMethodData().getNamespace();

  const std::string baseSocketName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_States_";

  std::string socketName = baseSocketName + "DATA";
  DataHandle<RealVector> interface = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();
  socketName = baseSocketName + "ISACCEPTED";
  DataHandle<CFreal> isAccepted = _sockets.getSocketSink<CFreal>(socketName)->getDataHandle();
  socketName = baseSocketName + "COORD";
  DataHandle<RealVector> coordinates = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();

  ///Check that the datahandle has same size as the TRS
  cf_assert(isAccepted.size() == trsStates->size());

  CFuint idx = 0;
  CFuint i = 0;
  for (CFuint iState = 0; iState < trsStates->size(); ++iState) {
    const CFuint stateID = (*trsStates)[iState];
    if(isAccepted[i] >= 0.){
      if (!isUpdated[stateID])
      {
//          CFout << "Updating: " << states[stateID]->getCoordinates() << "\n";
//          CFout << "...with data["<<idx<<"]: " << interface[idx] << "\n";
//          CFout << "...corresponding with: " << coordinates[idx] << "\n";

        *(states[stateID]) = interface[idx];

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) = 0.;
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
      }
      idx++;
    }
    else
    {
      if (!isUpdated[stateID])
      {
        _vFunction.evaluate(states[stateID]->getCoordinates(),*(states[stateID]));

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) = 0.;
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
      }
    }
    i++;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

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
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];
    const std::string baseSocketName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_States_";

    std::string socketName = baseSocketName + "DATA";
    _sockets.createSocketSink<RealVector>(socketName);
    socketName = baseSocketName + "ISACCEPTED";
    _sockets.createSocketSink<CFreal>(socketName);
    socketName = baseSocketName + "COORD";
    _sockets.createSocketSink<RealVector>(socketName);
  }

}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD
