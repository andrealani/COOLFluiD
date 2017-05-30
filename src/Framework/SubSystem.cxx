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
#include "Common/PE.hh"
#include "Common/CFPrintContainer.hh"
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
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void SubSystem::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("CreateNullMethods","Methods that where not configured are created Null.");
  options.addConfigOption< std::vector<string> >("Namespaces","Namespaces to be created in this SubSystem.");
  options.addConfigOption< std::vector<string> >
  ("Ranks","MPI ranks for each namespace (starting from 0) in the form START0:END0 START1:END1 etc.");
}

//////////////////////////////////////////////////////////////////////////////

SubSystem::SubSystem(const string& name)
  : ConfigObject(name),
    m_ranksCounter(0),
    m_namespaces(),
    m_has_null_methods(true),
    m_nbmethods(0)
{
  addConfigOptionsTo(this);

  setParameter("CreateNullMethods",&m_has_null_methods);
  setParameter("Namespaces",&m_namespaces);
  setParameter("Ranks",&m_ranks);
  
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
  
  // typedef std::vector<Common::SafePtr<Namespace> > NspVec;
  // NamespaceSwitcher& nsw = NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName());
  // NspVec lst = nsw.getAllNamespaces();
  // cf_assert(!lst.empty()); // there should be some models
  // const int rank = Common::PE::GetPE().GetRank("Default");
  // SafePtr<MeshData> currMD = CFNULL;
  // string mdName = "";
  // for(NspVec::iterator nsp = lst.begin(); nsp != lst.end(); ++nsp)
  // {
  //   const string nspaceName = (*nsp)->getName();
  //   if (PE::GetPE().isRankInGroup(rank, nspaceName)) {
  //    mdName = (*nsp)->getMeshDataName();
  //    break;
  //   }
  // }

  // vector <Common::SafePtr<MeshData> > meshDataVector = MeshDataStack::getInstance().getAllEntries();
  // for (CFuint iMeshData = 0; iMeshData < meshDataVector.size(); iMeshData++) {
  //   if (meshDataVector[iMeshData]->getName() == mdName) {
  //     MeshDataStack::getInstance().deleteEntry(mdName);
  //     break; 
  //   } 
  // }

  // delete all the SubSystemStatus
  SubSystemStatusStack::getInstance().deleteAllEntries();

  // delete all the namespaces
  NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).deleteAllNamespaces();

  deletePtr( m_int_param_reader );
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::configureNamespaces(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

  // configure just one namespace if none has been explicitly defined
  if (m_namespaces.empty())  
  {
    const string name = "Default";
    m_namespaces.push_back(name);
    Common::SafePtr<Namespace> ptr = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).createUniqueNamespace(name);
    configureNested(*ptr,args);
  
    // -----------------------------------------------------------//
    // AL: in order to enable multiple groups, each subsystem     //
    //     must define namespaces with unique names and different //
    //     from "Default" explicitly in the CFcase file           //
    // -----------------------------------------------------------//
    
    // The following can be enabled only if <subsystem name + namespace> 
    // is used as key for PE::GetPE() functions instead of namespace 
    // (e.g. "SubSyA.Default") 
    
    // PE::createGroup(name, name, granks, true);
    
    // make a unique list of all the ranks used in this subsystem
    vector<int> subSystemRanksUnique;
    const CFuint nbProc = PE::GetPE().GetProcessorCount("Default");  
    subSystemRanksUnique.resize(nbProc);
    for (CFuint r = 0; r < nbProc; ++r) {
      subSystemRanksUnique[r] = r;
    }
    
    PE::GetPE().createGroup("Default", SubSystemStatusStack::getCurrentName(), subSystemRanksUnique, true); 
  }
  // or configure the whole list that was defined
  else {
    // initial screening for namespaces to check for "|"
    // this is a hack which should be fixed in a cleaner way
    // namespaces can be already multiplied by using the "@" key
    // the problem is for the ranks ...
    vector<string> nsp;
    vector<string> ranks;
    vector<string> nspMulti;
    for (CFuint i = 0; i < m_namespaces.size(); ++i) {
      const size_t found = m_namespaces[i].find("|");
      if (found != string::npos) {
	nspMulti.push_back(m_namespaces[i]);
	vector<string> nsp2N =  StringOps::getWords(m_namespaces[i], '|');
	const string rootNsp = nsp2N[0];
	const CFuint nbNsp = StringOps::from_str<CFuint>(nsp2N[1]);
	vector<string> startEndRank = StringOps::getWords(m_ranks[i], ':');
	CFuint start = StringOps::from_str<CFuint>(startEndRank[0]);
	const CFuint end = StringOps::from_str<CFuint>(startEndRank[1]);
	const CFuint nbRanksPerNsp = nbNsp/(end+1-start);
	for (CFuint n = 0; n < nbNsp; ++n) {
	  const string nspPlusN = rootNsp + StringOps::to_str(n);
	  nsp.push_back(nspPlusN);
	  const string rankN = StringOps::to_str(start) + ":" + StringOps::to_str(start+nbRanksPerNsp-1);
	  ranks.push_back(rankN);
	  start += nbRanksPerNsp;
	}
      }
      else {
	nsp.push_back(m_namespaces[i]);
	// AL: this gory fix needs more testing
	if (m_ranks.size() > 0) {
	  cf_assert(i < m_ranks.size());
	  ranks.push_back(m_ranks[i]);
	}
      }
    }
    
    if (nsp.size() > m_namespaces.size()) {
      m_namespaces.resize(nsp.size());
      m_namespaces = nsp;
      CFLog(INFO, CFPrintContainer<vector<string> >("Namespaces = ", &m_namespaces) << "\n");
      m_ranks.resize(ranks.size());
      m_ranks = ranks;
      CFLog(INFO, CFPrintContainer<vector<string> >("Ranks = ", &m_ranks) << "\n");
    }
    
    CFuint counter = 0; 
    CFuint nbRanks = 0;
    std::vector<string>::iterator itr = m_namespaces.begin();
    vector<int> subSystemRanks;
    for(; itr != m_namespaces.end(); ++itr, ++counter) {
      Common::SafePtr<Namespace> ptr = NamespaceSwitcher::getInstance
	(SubSystemStatusStack::getCurrentName()).createUniqueNamespace(*itr);
      configureNested(*ptr,args);
      
      if (m_ranks.size() > 0 && (*itr != "Default")) {
	cf_assert(m_ranks.size() == m_namespaces.size());
	
	// -----------------------------------------------------------//
	// AL: in order to enable multiple groups, each subsystem     //
	//     must define namespaces with unique names and different //
	//     from "Default" explicitly in the CFcase file           //
	// -----------------------------------------------------------//
	
	// The following can be enabled only if <subsystem name + namespace> 
	// is used as key for PE::GetPE() functions instead of namespace 
	// (e.g. "SubSyA.Default") 
	
	// create corresponding MPI group
	vector<int> granks;
	fillGroupRanks(m_ranks[counter], granks);
	CFLog(VERBOSE, "SubSystem::configureNamespaces() => before creating group " << *itr << "\n");
	const string msg = "Ranks for group [" +  *itr + "] = ";
	CFLog(VERBOSE, CFPrintContainer<vector<int> >(msg, &granks));
	PE::GetPE().createGroup(*itr, *itr, granks, true);
	CFLog(VERBOSE, "SubSystem::configureNamespaces() => after  creating group " << *itr << "\n");
	copy(granks.begin(), granks.end(), back_inserter(subSystemRanks));
	nbRanks += granks.size();
      }
    }
    
    
    // make a unique list of all the ranks used in this subsystem
    vector<int> subSystemRanksUnique;
    
    if (nbRanks > 0) {
      sort(subSystemRanks.begin(), subSystemRanks.end());
      unique_copy(subSystemRanks.begin(), subSystemRanks.end(), back_inserter(subSystemRanksUnique));
      cf_assert(subSystemRanksUnique.size() <= PE::GetPE().GetProcessorCount("Default"));
    }
    
    // build a group, subset of MPI_COMM_WORLD corresponding to all ranks involved in this subsystem
    if (nbRanks == 0) {
      const CFuint nbProc = PE::GetPE().GetProcessorCount("Default");  
      subSystemRanksUnique.resize(nbProc);
      for (CFuint r = 0; r < nbProc; ++r) {
        subSystemRanksUnique[r] = r;
      }
    }
    
    PE::GetPE().createGroup("Default", SubSystemStatusStack::getCurrentName(), subSystemRanksUnique, true); 
  }
  
  if (m_ranksCounter == 0) {m_ranksCounter = PE::GetPE().GetProcessorCount("Default");}
}
    
//////////////////////////////////////////////////////////////////////////////

void SubSystem::configureSingletons( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  typedef std::vector<Common::SafePtr<Namespace> > NspVec;
  NspVec lst = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getAllNamespaces();

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
  NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).pushNamespace((*lst.begin())->getName());
  setEnableNamespaces(false);
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::setEnableNamespaces(const bool enable)
{
  NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).setEnabled(enable);
  MeshDataStack::getInstance().setEnabled(enable);
  PhysicalModelStack::getInstance().setEnabled(enable);
  SubSystemStatusStack::getInstance().setEnabled(enable);
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::configureNamespaceSingletons(Config::ConfigArgs& args, Common::SafePtr<Namespace> nsp)
{
  CFAUTOTRACE;

  cf_assert(nsp.isNotNull());

  const string meshDataName      = nsp->getMeshDataName();
  const string physicalModelName = nsp->getPhysicalModelName();
  const string physicalModelType = nsp->getPhysicalModelType();
  const string sysStatusName     = nsp->getSubSystemStatusName();

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
    // md->reallocate();
    md->setFactoryRegistry(getFactoryRegistry());
    configureNested(*md,args);
  }
  
  cf_assert(!physicalModelName.empty());
  Common::SafePtr<PhysicalModel> pm = PhysicalModelStack::getInstance().createUnique(physicalModelName);
  pm->setFactoryRegistry(getFactoryRegistry());
  if (!pm->isConfigured()) {
    configureNested(*pm,args);
  }

  cf_assert(!sysStatusName.empty());
  Common::SafePtr<SubSystemStatus> ss = SubSystemStatusStack::getInstance().createUnique(sysStatusName);
  ss->setFactoryRegistry(getFactoryRegistry());
  if (ss->getSubSystemName().empty()) {
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

  configureNamespaces(


args);
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
  
  vector<Method*> allMethods = getMethodList();
  vector<Method*>::iterator mitr = allMethods.begin();
  for ( ; mitr != allMethods.end(); ++mitr) {
    (*mitr)->setCommandGroups();
  }

  vector <Common::SafePtr<MeshData> > meshDataVector = 
    MeshDataStack::getInstance().getAllEntries();
  
  // loop on all mesh datas
  const int rank = PE::GetPE().GetRank("Default");
  vector <Common::SafePtr<MeshData> >::iterator mdata_itr = meshDataVector.begin();
  for( ; mdata_itr != meshDataVector.end(); ++mdata_itr)
  {
    Common::SafePtr<MeshData> mdata = *mdata_itr;
    if (PE::GetPE().isRankInGroup(rank, mdata->getPrimaryNamespace())) {
      vector< SafePtr<TopologicalRegionSet> > meshDataTRS = mdata->getTrsList();
      // all commands in the namespace of this mesh data
      vector< SafePtr<NumericalCommand> > meshDataNamespaceCommands;
      
      // loop on all methods
      vector<Method*>::iterator mitr = allMethods.begin();
      for ( ; mitr != allMethods.end(); ++mitr) {
	if( mdata->match ( (*mitr)->getNamespace() ) )
	  putCommandsFromMethod( *mitr, meshDataNamespaceCommands );
      }
      
      // in this MeshData set the TRSs pointers in the commands
      CommandsToTRSMapper mapper ( meshDataNamespaceCommands, meshDataTRS );
      mapper.mapComsToTrs();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubSystem::registActionListeners()
{
  CFAUTOTRACE;
  
  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  
  const string ssname = SubSystemStatusStack::getCurrentName();   
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_PLUGSOCKETS"),
			     this,&SubSystem::allocateAllSocketsAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_BUILDPHYSICALMODEL"), 
			     this,&SubSystem::buildPhysicalModelAction);
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

Common::Signal::return_t SubSystem::buildPhysicalModelAction(Common::Signal::arg_t eBuild)
{
  buildPhysicalModel();
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

void SubSystem::fillGroupRanks(const string rankString, vector<int>& granks)
{
  // ranks are in the form START:END
  vector<string> ranks = StringOps::getWords(rankString, ':');
  cf_assert(ranks.size() == 2);
  const CFuint start = StringOps::from_str<CFuint>(ranks[0]);
  const CFuint end   = StringOps::from_str<CFuint>(ranks[1]);
  const CFuint gsize = end-start+1; 
  m_ranksCounter += end-start+1;
  cf_assert(m_ranksCounter > 0);
  
  granks.resize(gsize);
  for (CFuint r = 0; r < gsize; ++r) {
    granks[r] = start+r;
  }
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
