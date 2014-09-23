#include "Common/EventHandler.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/SubSystemStatus.hh"

#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"
#include "MeshAdapterSpringAnalogy/SpringAnalogy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  using namespace Common;
  using namespace Framework;
  using namespace Common;

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SpringAnalogy,
               MeshAdapterMethod,
               MeshAdapterSpringAnalogyModule,
               1>
springAnalogyMeshAdapterMethodProvider("SpringAnalogy");

//////////////////////////////////////////////////////////////////////////////

void SpringAnalogy::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareComds","Types of the initializing commands.");
   options.addConfigOption< std::string >("UpdateMeshCom","Update Mesh command to run.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::vector<std::string> >("PrepareNames","Names of the initializing commands.");
}

//////////////////////////////////////////////////////////////////////////////

SpringAnalogy::SpringAnalogy(const std::string& name)
  : MeshAdapterMethod(name),
  _prepareTypeStr(0),
  _prepareNameStr(0)
{
   addConfigOptionsTo(this);

  _data.reset(new SpringAnalogyData(this));

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

SpringAnalogy::~SpringAnalogy()
{
  clearPrepareComs();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> SpringAnalogy::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void SpringAnalogy::clearPrepareComs()
{
  _prepares.clear();
  vector<SelfRegistPtr<SpringAnalogyCom> >().swap(_prepares);
}


//////////////////////////////////////////////////////////////////////////////

void SpringAnalogy::configure ( Config::ConfigArgs& args )
{
  MeshAdapterMethod::configure(args);
  configureNested ( _data.getPtr(), args );

  // add configures to the SpringAnalogyCom's

  configureCommand<SpringAnalogyData,SpringAnalogyComProvider>( args, _setup,_setupStr,_data);

  configureCommand<SpringAnalogyData,SpringAnalogyComProvider>( args, _unSetup,_unSetupStr,_data);

  configureCommand<SpringAnalogyData,SpringAnalogyComProvider>( args, _transformMesh,_transformMeshStr,_data);

//  configureCommand<SpringAnalogyData,SpringAnalogyComProvider>( args, _prepare,_prepareStr,_data);

  cf_assert(_prepareTypeStr.size() == _prepareNameStr.size());

  _prepares.resize(_prepareTypeStr.size());
  for(CFuint i = 0; i < _prepares.size(); ++i) {

    CFLog(VERBOSE, "PREPARE type = " << _prepareTypeStr[i] << CFendl);
    CFLog(VERBOSE, "PREPARE name = " << _prepareNameStr[i] << CFendl);

    configureCommand<SpringAnalogyCom,
                     SpringAnalogyData,
                     SpringAnalogyComProvider>( args, _prepares[i],
                                           _prepareTypeStr[i],
                                           _prepareNameStr[i],
                                           _data);

    cf_assert(_prepares[i].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpringAnalogy::setMethodImpl()
{
  MeshAdapterMethod::setMethodImpl();

  _data->getGeoWithNodesBuilder()->setup();

  setupCommandsAndStrategies();
  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SpringAnalogy::unsetMethodImpl()
{
  _unSetup->execute();
  unsetupCommandsAndStrategies();

  MeshAdapterMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void SpringAnalogy::adaptMeshImpl()
{
  CFAUTOTRACE;

  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  std::string msg;

  if(!(SubSystemStatusStack::getActive()->getNbIter() % _adaptRate))
  {
    (_data.getPtr())->resetCurrentStep();
    const CFuint nbSteps = (_data.getPtr())->getTotalNbSteps();

    // raise the event for the backup of data prior to change
    event_handler->call_signal ( "CF_ON_MESHADAPTER_BEFOREMESHUPDATE", msg );

    for(CFuint iStep = 1;iStep<nbSteps+1;iStep++)
    {
      std::cout << "********************** Step " << iStep<< std::endl;
      for(CFuint i=0;i<_prepares.size();i++)
      {
        cf_assert(_prepares[i].isNotNull());
        _prepares[i]->execute();
      }

      _transformMesh->execute();
      (_data.getPtr())->updateCurrentStep();
    }

    // Raise the event for the update of the data in the other methods
    event_handler->call_signal ( "CF_ON_MESHADAPTER_AFTERMESHUPDATE", msg );

  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

