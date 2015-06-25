// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/EventHandler.hh"
#include "Common/PE.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "SimpleGlobalMeshAdapter/SimpleMeshAdapter.hh"
#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  using namespace Common;
  using namespace Framework;
  using namespace Common;

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SimpleMeshAdapter,
               MeshAdapterMethod,
               SimpleGlobalMeshAdapterModule,
               1>
simpleMeshAdapterMeshAdapterMethodProvider("SimpleMeshAdapter");

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("RemeshCondition","Condition to remesh.");
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
   options.addConfigOption< std::string >("PrepareCom","Command to perform the preparatory steps before mesh adaptation.");
   options.addConfigOption< std::string >("MeshGeneratorCom","Command to generate the new mesh.");
   options.addConfigOption< std::string >("MeshReaderCom","Command to read the new mesh.");
   options.addConfigOption< std::string >("MeshInterpolatorCom","Command to interpolate the solution on the new mesh.");
   options.addConfigOption< std::string >("MeshWriterCom","Command to write the new mesh with the solution.");
   options.addConfigOption< std::string >("UnSetupCom","UnSetupCommand to run. This command seldomly needs overriding.");
}

//////////////////////////////////////////////////////////////////////////////

SimpleMeshAdapter::SimpleMeshAdapter(const std::string& name)
  : MeshAdapterMethod(name)
{
   addConfigOptionsTo(this);

  _data.reset(new SimpleMeshAdapterData(this));

  _setupStr = "StdSetup";
   setParameter("SetupCom",&_setupStr);

  _prepareStr = "StdPrepare";
   setParameter("PrepareCom", &_prepareStr);

  _runMeshGeneratorStr = "DummyMeshGenerator";
   setParameter("MeshGeneratorCom",&_runMeshGeneratorStr);

  _readNewMeshStr = "StdMeshReader";
   setParameter("MeshReaderCom",&_readNewMeshStr);

  _interpolateSolStr = "DummyMeshInterpolator";
   setParameter("MeshInterpolatorCom",&_interpolateSolStr);

  _writeInterpolatedMeshStr = "StdMeshWriter";
   setParameter("MeshWriterCom",&_writeInterpolatedMeshStr);

  _unSetupStr = "StdUnSetup";
   setParameter("UnSetupCom",&_unSetupStr);

  _remeshConditionStr = "NullRemeshCondition";
   setParameter("RemeshCondition",&_remeshConditionStr);
}

//////////////////////////////////////////////////////////////////////////////

SimpleMeshAdapter::~SimpleMeshAdapter()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> SimpleMeshAdapter::getMethodData() const
{
  return _data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapter::configure ( Config::ConfigArgs& args )
{

  MeshAdapterMethod::configure(args);
  configureNested ( _data.getPtr(), args );

  // add configures to the SimpleMeshAdapterCom's
  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _setup,_setupStr,_data);

  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _prepare,_prepareStr,_data);

  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _runMeshGenerator,_runMeshGeneratorStr,_data);

  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _readNewMesh,_readNewMeshStr,_data);

  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _interpolateSol, _interpolateSolStr, _data);

  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _writeInterpolatedMesh, _writeInterpolatedMeshStr, _data);

  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _unSetup,_unSetupStr,_data);

  configureCommand<SimpleMeshAdapterData,SimpleMeshAdapterComProvider>( args, _remeshCondition,_remeshConditionStr
,_data);

}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapter::setNonRootMethods()
{
  MultiMethodHandle<MeshCreator> meshCreator = getMeshCreator();
  cf_assert(meshCreator.size() == 1);

  MultiMethodHandle<OutputFormatter> outputFormat = getOutputFormatter();
  cf_assert(outputFormat.size() == 1);

  (meshCreator[0])->setIsNonRootMethod();
  (outputFormat[0])->setIsNonRootMethod();
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapter::setMethodImpl()
{
  MeshAdapterMethod::setMethodImpl();

  _data->setSpaceMethod(getSpaceMethod());
  _data->setConvergenceMethod(getConvergenceMethod());
  _data->setMeshCreator(getMeshCreator());
  _data->setOutputFormatter(getOutputFormatter());

  setupCommandsAndStrategies();
  _setup->execute();
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapter::unsetMethodImpl()
{
  _unSetup->execute();
  unsetupCommandsAndStrategies();

  MeshAdapterMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapter::adaptMeshImpl()
{
  //do nothing here
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapter::remeshImpl()
{

  CFAUTOTRACE;

  const std::string otherNameSpace = _data->getOtherNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(otherNameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
  
  _remeshCondition->execute();
  CFout <<" Need remeshing ?? : "<< _data->isNeedRemeshing() <<"\n";

  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  std::string msg;
  const std::string ssname = SubSystemStatusStack::getCurrentName();   
  
  if(!(subsystemStatus->getNbIter() % _adaptRate) && _data->isNeedRemeshing())
  {
    // Raise the event for the backup of data prior to change
    //    event_handler->call_signal (event_handler->key(ssname,"CF_ON_MESHADAPTER_BEFOREGLOBALREMESHING"), msg );

    _prepare->execute();
    _runMeshGenerator->execute();
    _readNewMesh->execute();
    _interpolateSol->execute();
    _writeInterpolatedMesh->execute();
    
   //   PE::GetPE().setBarrier();
   //   PE::GetPE().setBarrier();
    
    //Raise the event for the update of the data in the other methods
    event_handler->call_signal (event_handler->key(ssname, "CF_ON_MESHADAPTER_AFTERGLOBALREMESHING"), msg );
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

