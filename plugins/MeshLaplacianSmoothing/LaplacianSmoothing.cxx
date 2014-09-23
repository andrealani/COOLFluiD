#include "Common/EventHandler.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/SubSystemStatus.hh"

#include "MeshLaplacianSmoothing/MeshLaplacianSmoothing.hh"
#include "MeshLaplacianSmoothing/LaplacianSmoothing.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  using namespace Common;
  using namespace Framework;
  using namespace Common;

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LaplacianSmoothing,
               MeshAdapterMethod,
               MeshLaplacianSmoothingModule,
               1>
springAnalogyMeshAdapterMethodProvider("LaplacianSmoothing");

//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothing::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareComds","Types of the initializing commands.");
   options.addConfigOption< std::string >("UpdateMeshCom","Update Mesh command to run.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareNames","Names of the initializing commands.");
}

//////////////////////////////////////////////////////////////////////////////

LaplacianSmoothing::LaplacianSmoothing(const std::string& name)
  : MeshAdapterMethod(name),
  _prepareTypeStr(0),
  _prepareNameStr(0)
{
   addConfigOptionsTo(this);

  _data.reset(new LaplacianSmoothingData(this));

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

LaplacianSmoothing::~LaplacianSmoothing()
{
  clearPrepareComs();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> LaplacianSmoothing::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothing::clearPrepareComs()
{
  _prepares.clear();
  vector<SelfRegistPtr<LaplacianSmoothingCom> >().swap(_prepares);
}


//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothing::configure ( Config::ConfigArgs& args )
{
  MeshAdapterMethod::configure(args);
  configureNested ( _data.getPtr(), args );

  // add configures to the LaplacianSmoothingCom's

  configureCommand<LaplacianSmoothingData,LaplacianSmoothingComProvider>( args, _setup,_setupStr,_data);

  configureCommand<LaplacianSmoothingData,LaplacianSmoothingComProvider>( args, _unSetup,_unSetupStr,_data);

  configureCommand<LaplacianSmoothingData,LaplacianSmoothingComProvider>( args, _transformMesh,_transformMeshStr,_data);

//  configureCommand<LaplacianSmoothingData,LaplacianSmoothingComProvider>( args, _prepare,_prepareStr,_data);

  cf_assert(_prepareTypeStr.size() == _prepareNameStr.size());

  _prepares.resize(_prepareTypeStr.size());
  for(CFuint i = 0; i < _prepares.size(); ++i) {

    CFLog(VERBOSE, "PREPARE type = " << _prepareTypeStr[i] << CFendl);
    CFLog(VERBOSE, "PREPARE name = " << _prepareNameStr[i] << CFendl);

    configureCommand<LaplacianSmoothingCom,
                     LaplacianSmoothingData,
                     LaplacianSmoothingComProvider>( args, _prepares[i],
                                           _prepareTypeStr[i],
                                           _prepareNameStr[i],
                                           _data);

    cf_assert(_prepares[i].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothing::setMethodImpl()
{
  MeshAdapterMethod::setMethodImpl();

  setupCommandsAndStrategies();
  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothing::unsetMethodImpl()
{
  _unSetup->execute();
  unsetupCommandsAndStrategies();

  MeshAdapterMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void LaplacianSmoothing::adaptMeshImpl()
{

  CFAUTOTRACE;

  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();

  if(!(SubSystemStatusStack::getActive()->getNbIter() % _adaptRate))
  {
    std::string msg;

    //Raise the event for the backup of data prior to change
    event_handler->call_signal ( "CF_ON_MESHADAPTER_BEFOREMESHUPDATE", msg );

    for(CFuint i=0;i<_prepares.size();i++)
    {
      cf_assert(_prepares[i].isNotNull());
      _prepares[i]->execute();
    }

    _transformMesh->execute();

    // Raise the event for the update of the data in the other methods
    event_handler->call_signal ( "CF_ON_MESHADAPTER_AFTERMESHUPDATE", msg );
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

