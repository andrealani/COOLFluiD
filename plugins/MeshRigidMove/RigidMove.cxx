#include "Common/EventHandler.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/SubSystemStatus.hh"

#include "MeshRigidMove/RigidMove.hh"
#include "MeshRigidMove/MeshRigidMove.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  using namespace Common;
  using namespace Framework;
  using namespace Common;

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RigidMove,
               MeshAdapterMethod,
               MeshRigidMoveModule,
               1>
rigidMoveMeshAdapterMethodProvider("RigidMove");

//////////////////////////////////////////////////////////////////////////////

void RigidMove::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareComds","Types of the initializing commands.");
   options.addConfigOption< std::string >("UpdateMeshCom","Update Mesh command to run.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareNames","Names of the initializing commands.");
}

//////////////////////////////////////////////////////////////////////////////

RigidMove::RigidMove(const std::string& name)
  : MeshAdapterMethod(name),
  _prepareTypeStr(0),
  _prepareNameStr(0)
{
   addConfigOptionsTo(this);

  _data.reset(new RigidMoveData(this));

  _setupStr = "StdSetup";
  setParameter("SetupCom",&_setupStr);
  
  _unSetupStr = "StdUnSetup";
  setParameter("UnSetupCom",&_unSetupStr);
  
  _transformMeshStr = "UpdateMesh";
  setParameter("UpdateMeshCom",&_transformMeshStr);
  
  _prepareTypeStr.clear();
  setParameter("PrepareComds",&_prepareTypeStr);
  
  _prepareNameStr.clear();
  setParameter("PrepareNames",&_prepareNameStr);
}

//////////////////////////////////////////////////////////////////////////////

RigidMove::~RigidMove()
{
  CFLog(VERBOSE, "RigidMove::~RigidMove() START\n");
  clearPrepareComs();
  CFLog(VERBOSE, "RigidMove::~RigidMove() END\n");
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> RigidMove::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void RigidMove::clearPrepareComs()
{
  CFLog(VERBOSE, "RigidMove::clearPrepareComs() START\n");
  vector<SelfRegistPtr<RigidMoveCom> >().swap(_prepares);
  CFLog(VERBOSE, "RigidMove::clearPrepareComs() END\n");
}


//////////////////////////////////////////////////////////////////////////////

void RigidMove::configure ( Config::ConfigArgs& args )
{
  MeshAdapterMethod::configure(args); 
  configureNested ( _data.getPtr(), args );
  
  // add configures to the RigidMoveCom's
  configureCommand<RigidMoveData,RigidMoveComProvider>( args, _setup,_setupStr,_data);

  configureCommand<RigidMoveData,RigidMoveComProvider>( args, _unSetup,_unSetupStr,_data);

  configureCommand<RigidMoveData,RigidMoveComProvider>( args, _transformMesh,_transformMeshStr,_data);

  cf_assert(_prepareTypeStr.size() == _prepareNameStr.size());
  if(_prepareTypeStr.size() == 0)
  {
    _prepareTypeStr.resize(1);
    _prepareNameStr.resize(1);
    _prepareTypeStr[0]="StdPrepare";
    _prepareNameStr[0]="StdPrepare1";
  }

  _prepares.resize(_prepareTypeStr.size());

  for(CFuint i = 0; i < _prepares.size(); ++i) {

    CFLog(INFO, "PREPARE type = " << _prepareTypeStr[i] << "\n");
    CFLog(INFO, "PREPARE name = " << _prepareNameStr[i] << "\n");

    configureCommand<RigidMoveCom,
                     RigidMoveData,
                     RigidMoveComProvider>( args, _prepares[i],
                                           _prepareTypeStr[i],
                                           _prepareNameStr[i],
                                           _data);

    cf_assert(_prepares[i].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void RigidMove::setMethodImpl()
{
  MeshAdapterMethod::setMethodImpl();

  setupCommandsAndStrategies();
  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void RigidMove::unsetMethodImpl()
{
  CFLog(VERBOSE, "RigidMove::unsetMethodImpl() START\n");
  _unSetup->execute(); 
  unsetupCommandsAndStrategies();
  MeshAdapterMethod::unsetMethodImpl();
  CFLog(VERBOSE, "RigidMove::unsetMethodImpl() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void RigidMove::adaptMeshImpl()
{
  CFAUTOTRACE;

  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  std::string msg;
  const std::string ssname = SubSystemStatusStack::getCurrentName();   
  if(!(SubSystemStatusStack::getActive()->getNbIter() % _adaptRate))
  {
    // raise the event for the backup of data prior to change
    event_handler->call_signal (event_handler->key(ssname, "CF_ON_MESHADAPTER_BEFOREMESHUPDATE"), msg );

    for(CFuint i=0;i<_prepares.size();i++)
    {
      cf_assert(_prepares[i].isNotNull());
      _prepares[i]->execute();
    }
    _transformMesh->execute();

    // raise the event for the update of the data in the other methods
    event_handler->call_signal (event_handler->key(ssname,"CF_ON_MESHADAPTER_AFTERMESHUPDATE"), msg );
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RigidMove

  } // namespace Numerics

} // namespace COOLFluiD

