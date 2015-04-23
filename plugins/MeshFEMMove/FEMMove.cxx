#include "Common/EventHandler.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/SubSystemStatus.hh"

#include "MeshFEMMove/MeshFEMMove.hh"
#include "MeshFEMMove/FEMMove.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  using namespace Common;
  using namespace Framework;
  using namespace Common;

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<FEMMove,
               MeshAdapterMethod,
               MeshFEMMoveModule,
               1>
FEMMoveMeshAdapterMethodProvider("FEMMove");

//////////////////////////////////////////////////////////////////////////////

void FEMMove::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareComds","Types of the initializing commands.");
   options.addConfigOption< std::string >("UpdateMeshCom","Update Mesh command to run.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareNames","Names of the initializing commands.");
}

//////////////////////////////////////////////////////////////////////////////

FEMMove::FEMMove(const std::string& name)
  : MeshAdapterMethod(name),
  _prepareTypeStr(0),
  _prepareNameStr(0)
{
   addConfigOptionsTo(this);

  _data.reset(new FEMMoveData(this));

  _setupStr = "StdSetup";
   setParameter("SetupCom",&_setupStr);



  _unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&_unSetupStr);



  _transformMeshStr = "UpdateMesh";
   setParameter("UpdateMeshCom",&_transformMeshStr);



//   _prepareStr = "StdPrepare";
//   addConfigOption("PrepareCom",
//       "Command to prepare the solution before the iteration process.",
//       &_prepareStr);

  _prepareTypeStr.clear();
   setParameter("PrepareComds",&_prepareTypeStr);



  _prepareNameStr.clear();
   setParameter("PrepareNames",&_prepareNameStr);



}

//////////////////////////////////////////////////////////////////////////////

FEMMove::~FEMMove()
{
  clearPrepareComs();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> FEMMove::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void FEMMove::clearPrepareComs()
{
  _prepares.clear();
  vector<SelfRegistPtr<FEMMoveCom> >().swap(_prepares);
}


//////////////////////////////////////////////////////////////////////////////

void FEMMove::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  MeshAdapterMethod::configure(args);
  configureNested ( _data.getPtr(), args );

  // add configures to the FEMMoveCom's
  configureCommand<FEMMoveData,FEMMoveComProvider>( args, _setup,_setupStr,_data);

  configureCommand<FEMMoveData,FEMMoveComProvider>( args, _unSetup,_unSetupStr,_data);

  configureCommand<FEMMoveData,FEMMoveComProvider>( args, _transformMesh,_transformMeshStr,_data);

//  configureCommand<FEMMoveData,FEMMoveComProvider>( args, _prepare,_prepareStr,_data);

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

    configureCommand<FEMMoveCom,
                     FEMMoveData,
                     FEMMoveComProvider>( args, _prepares[i],
                                           _prepareTypeStr[i],
                                           _prepareNameStr[i],
                                           _data);

    cf_assert(_prepares[i].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void FEMMove::setNonRootMethods()
{
  MultiMethodHandle<ConvergenceMethod> convergenceMethod = getConvergenceMethod();
  cf_assert(convergenceMethod.size() == 1);

  (convergenceMethod[0])->setIsNonRootMethod();
}

//////////////////////////////////////////////////////////////////////////////

void FEMMove::setMethodImpl()
{
  MeshAdapterMethod::setMethodImpl();

  _data->setConvergenceMethod(getConvergenceMethod());

  setupCommandsAndStrategies();

  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void FEMMove::unsetMethodImpl()
{
  _unSetup->execute();

  unsetupCommandsAndStrategies();

  MeshAdapterMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void FEMMove::adaptMeshImpl()
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
    event_handler->call_signal (event_handler->key(ssname, "CF_ON_MESHADAPTER_AFTERMESHUPDATE"), msg );
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FEMMove

  } // namespace Numerics

} // namespace COOLFluiD

