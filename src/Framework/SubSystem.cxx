// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <functional>

#include "Common/COOLFluiD.hh"

#ifdef CF_HAVE_VALGRIND
  #include <valgrind/callgrind.h>
#else // no valgrind so define empty macros
  #define CALLGRIND_START_INSTRUMENTATION
  #define CALLGRIND_STOP_INSTRUMENTATION
#endif // CF_HAVE_VALGRIND

#include "Common/EventHandler.hh"
#include "Framework/SubSystem.hh"
#include "Framework/Method.hh"
#include "Framework/NumericalCommand.hh"
#include "Framework/NumericalStrategy.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/CommandsToTRSMapper.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/InteractiveParamReader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void SubSystem::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("CreateNullMethods","Methods that where not configured are created Null.");
   options.addConfigOption< std::vector<std::string> >("Namespaces","Namespaces to be created in this SubSystem.");
}

//////////////////////////////////////////////////////////////////////////////

SubSystem::SubSystem(const std::string& name)
  : ConfigObject(name),
    m_namespaces(),
    m_has_null_methods(true),
    m_nbmethods(0)
{
  addConfigOptionsTo(this);

  setParameter("CreateNullMethods",&m_has_null_methods);
  setParameter("Namespaces",&m_namespaces);

  m_int_param_reader = new InteractiveParamReader("InteractiveParamReader");
}

//////////////////////////////////////////////////////////////////////////////

SubSystem::~SubSystem()
{
  CFAUTOTRACE;

  // delete all the PhysicalModels
  PhysicalModelStack::getInstance().deleteAllEntries();

  // delete all the MeshDatas
  MeshDataStack::getInstance().deleteAllEntries();

  // delete all the SubSystemStatus
  SubSystemStatusStack::getInstance().deleteAllEntries();

  // delete all the namespaces
  NamespaceSwitcher::getInstance().deleteAllNamespaces();

  deletePtr( m_int_param_reader );
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::configureNamespaces(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

  // configure just one namespace if none has been explicitly defined
  if (m_namespaces.empty())
  {
    std::string name = "Default";
    m_namespaces.push_back(name);
    Common::SafePtr<Namespace> ptr = NamespaceSwitcher::getInstance().createUniqueNamespace(name);
    configureNested(*ptr,args);
  }
  // or configure the whole list that was defined
  else
  {
    std::vector<std::string>::iterator itr = m_namespaces.begin();
    for(; itr != m_namespaces.end(); ++itr)
    {
      Common::SafePtr<Namespace> ptr = NamespaceSwitcher::getInstance().createUniqueNamespace(*itr);
      configureNested(*ptr,args);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::configureSingletons( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  typedef std::vector<Common::SafePtr<Namespace> > NspVec;
  NspVec lst = NamespaceSwitcher::getInstance().getAllNamespaces();

  cf_assert(!lst.empty());

  // process the first Namespace which is the Default one
  // if none was configured, create a Default Singleton of each type
  NspVec::iterator nsp = lst.begin();

  if ((*nsp)->getMeshDataName().empty()) (*nsp)->setMeshDataName("Default");
  if ((*nsp)->getPhysicalModelName().empty()) (*nsp)->setPhysicalModelName((*nsp)->getPhysicalModelType());
  if ((*nsp)->getSubSystemStatusName().empty()) (*nsp)->setSubSystemStatusName("SubSystemStatus");

  for(; nsp != lst.end(); ++nsp) {
    configureNamespaceSingletons(args,*nsp);
  }

  // Push the Default namespace into the stack
  setEnableNamespaces(true);
  NamespaceSwitcher::getInstance().pushNamespace((*lst.begin())->getName());
  setEnableNamespaces(false);
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::setEnableNamespaces(const bool enable)
{

  NamespaceSwitcher::getInstance().setEnabled(enable);
  MeshDataStack::getInstance().setEnabled(enable);
  PhysicalModelStack::getInstance().setEnabled(enable);
  SubSystemStatusStack::getInstance().setEnabled(enable);

}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::configureNamespaceSingletons(Config::ConfigArgs& args, Common::SafePtr<Namespace> nsp)
{
  CFAUTOTRACE;

  cf_assert(nsp.isNotNull());

  const std::string meshDataName      = nsp->getMeshDataName();
  const std::string physicalModelName = nsp->getPhysicalModelName();
  const std::string physicalModelType = nsp->getPhysicalModelType();
  const std::string sysStatusName     = nsp->getSubSystemStatusName();

  CFLog(NOTICE,"-------------------------------------------------------------\n");
  CFLog(NOTICE,"Setting Namespace : " << nsp->getName() << "\n");
  CFLog(NOTICE,"MeshData          : " << meshDataName << "\n");
  CFLog(NOTICE,"PhysicalModelName : " << physicalModelName << "\n");
  CFLog(NOTICE,"PhysicalModelType : " << physicalModelType << "\n");
  CFLog(NOTICE,"SubSysStatus      : " << sysStatusName << "\n");
  CFLog(NOTICE,"-------------------------------------------------------------\n");

  cf_assert(!meshDataName.empty());
  Common::SafePtr<MeshData> md = MeshDataStack::getInstance().createUnique(meshDataName);

  if (!md->isConfigured()) {
    md->reallocate();
    configureNested(*md,args);
  }

  cf_assert(!physicalModelName.empty());
  Common::SafePtr<PhysicalModel> pm = PhysicalModelStack::getInstance().createUnique(physicalModelName);

  if (!pm->isConfigured()) {
    configureNested(*pm,args);
  }

  cf_assert(!sysStatusName.empty());
  Common::SafePtr<SubSystemStatus> ss = SubSystemStatusStack::getInstance().createUnique(sysStatusName);

  if (ss->getSubSystemName().empty())
  {
    ss->setSubSystemName(getName());
    ss->setMovingMesh(false);
    configureNested(*ss,args);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ConfigObject::configure(args);

  configureNamespaces(args);
  configureSingletons(args);

  configureNested ( m_int_param_reader, args );
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::allocateAllSockets()
{
  CFAUTOTRACE;

  // when multiple MeshData's are implemented, loop over them
  typedef std::vector<Common::SafePtr<MeshData> > MDLst;
  MDLst mds = MeshDataStack::getInstance().getAllEntries();
  MDLst::iterator meshData = mds.begin();
  for(; meshData != mds.end(); ++meshData)
  {
     (*meshData)->allocateConnectivity();
//     (*meshData)->allocateSockets();
  }

#if 0
  this->allocateSockets();
#endif
  checkAllSockets();
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::deallocateAllSockets()
{
  CFAUTOTRACE;

#if 0
  this->deallocateSockets();
#endif

  // when multiple MeshData's are implemented, loop over them
  typedef std::vector<Common::SafePtr<MeshData> > MDLst;
  MDLst mds = MeshDataStack::getInstance().getAllEntries();
  MDLst::iterator meshData = mds.begin();
  for(; meshData != mds.end(); ++meshData)
  {
    (*meshData)->deallocateNodesStates();
    (*meshData)->deallocateConnectivity();
  }

}

//////////////////////////////////////////////////////////////////////////////

std::pair < SubSystem::NCommandVec_t, SubSystem::NStrategyVec_t > SubSystem::getAllCommandsAndStrategies() // const
{
  NCommandVec_t  allComs;
  NStrategyVec_t allStrats;

  vector<Method*> allMethods = getMethodList();

  vector<Method*>::iterator mitr = allMethods.begin();
  for ( ; mitr != allMethods.end(); ++mitr )
  {
    putCommandsFromMethod ( *mitr, allComs );
    putStrategiesFromMethod ( *mitr, allStrats );
  }
  return std::pair < NCommandVec_t, NStrategyVec_t > (allComs,allStrats);
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::putCommandsFromMethod(Method* const mtd, SubSystem::NCommandVec_t& coms) // const
{
  vector<Common::SafePtr<NumericalCommand> > someComs = mtd->getCommandList();
  unique_copy(someComs.begin(),someComs.end(),back_inserter(coms));
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::putStrategiesFromMethod(Method* const mtd, SubSystem::NStrategyVec_t& strat) // const
{
  vector<Common::SafePtr<NumericalStrategy> > someStrat = mtd->getStrategyList();
  unique_copy(someStrat.begin(),someStrat.end(),back_inserter(strat));
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::setCommands()
{
  CFAUTOTRACE;

//   std::pair< NCommandVec_t , NStrategyVec_t > allCommandsStrategies = getAllCommandsAndStrategies();

  vector<Method*> allMethods = getMethodList();

  vector<Method*>::iterator mitr = allMethods.begin();
  for ( ; mitr != allMethods.end(); ++mitr )
  {
    (*mitr)->setCommandGroups();
  }

  vector <Common::SafePtr<MeshData> > meshDataVector  = MeshDataStack::getInstance().getAllEntries();

  // loop on all mesh datas
  vector <Common::SafePtr<MeshData> >::iterator mdata_itr = meshDataVector.begin();
  for( ; mdata_itr != meshDataVector.end(); ++mdata_itr)
  {
    Common::SafePtr<MeshData> mdata = *mdata_itr;
    vector< SafePtr<TopologicalRegionSet> > meshDataTRS = mdata->getTrsList();

    vector< SafePtr<NumericalCommand> > meshDataNamespaceCommands; // all commands in the namespace of this mesh data

    // loop on all methods
    vector<Method*>::iterator mitr = allMethods.begin();
    for ( ; mitr != allMethods.end(); ++mitr )
    {
      if( mdata->match ( (*mitr)->getNamespace() ) )
        putCommandsFromMethod( *mitr, meshDataNamespaceCommands );
    }

    // in this MeshData set the TRSs pointers in the commands
    CommandsToTRSMapper mapper ( meshDataNamespaceCommands, meshDataTRS );
    mapper.mapComsToTrs();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::registActionListeners()
{
  CFAUTOTRACE;
  
  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  
  const std::string ssname = SubSystemStatusStack::getCurrentName();   
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_PLUGSOCKETS"),
			     this,&SubSystem::allocateAllSocketsAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_BUILDMESHDATA"), 
			     this,&SubSystem::buildMeshDataAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_SETUP"),
			     this,&SubSystem::setupAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_RUN"),
			     this,&SubSystem::runAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_UNPLUGSOCKETS"), 
			     this,&SubSystem::deallocateAllSocketsAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_UNSETUP"),
			     this,&SubSystem::unsetupAction);
}
    
//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SubSystem::allocateAllSocketsAction(Common::Signal::arg_t eSocketsPlug)
{
  allocateAllSockets();
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SubSystem::deallocateAllSocketsAction(Common::Signal::arg_t eSocketsUnPlug)
{
  deallocateAllSockets();
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SubSystem::buildMeshDataAction(Common::Signal::arg_t eBuild)
{
  buildMeshData();
  return Common::Signal::return_t ();
}
    
//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SubSystem::setupAction(Common::Signal::arg_t eModifyRestart)
{
  setup();
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SubSystem::runAction(Common::Signal::arg_t eModifyRestart)
{
  CALLGRIND_START_INSTRUMENTATION;
  run();
  CALLGRIND_STOP_INSTRUMENTATION;
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SubSystem::unsetupAction(Common::Signal::arg_t eModifyRestart)
{
  unsetup();
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::allocateSockets()
{
  CFAUTOTRACE;

  vector<Method*> methodList = getMethodList();
  for_each(methodList.begin(),
           methodList.end(),
           mem_fun(&Method::allocateMethodSockets));
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::deallocateSockets()
{
  CFAUTOTRACE;

  vector<Method*> methodList = getMethodList();
  for_each(methodList.begin(),
           methodList.end(),
           mem_fun(&Method::deallocateMethodSockets));
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::checkAllSockets()
{
  CFAUTOTRACE;

  vector<Method*> methodList = getMethodList();
  for_each(methodList.begin(),
           methodList.end(),
           mem_fun(&Method::checkAllSockets));
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::setParentNamespaceInMethodSockets()
{
  CFAUTOTRACE;

  vector<Method*> methodList = getMethodList();

  // sets the Parent Namespace in all sockets
  for_each(methodList.begin(),
           methodList.end(),
           mem_fun(&Method::setSocketNamespaces));
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
