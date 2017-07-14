// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Common/CFLog.hh"
#include "Common/MemFunArg.hh"
#include "Common/EventHandler.hh"
#include "Common/NullPointerException.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/CFEnv.hh"
#include "Environment/CFEnvVars.hh"

#include "Framework/OnlyMeshSubSystem.hh"
#include "Framework/MeshData.hh"
#include "Framework/CFL.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/ErrorEstimatorMethod.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/OutputFormatter.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/InteractiveParamReader.hh"
#include "Framework/CommandsToTRSMapper.hh"
#include "Framework/PathAppender.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/Namespace.hh"
#include "Framework/Framework.hh"
#include "Framework/SimulationStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::filesystem;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<OnlyMeshSubSystem,
               SubSystem,
               FrameworkLib,
               1>
onlyMeshSubSysProvider("OnlyMeshSubSystem");

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("MeshAdapterNames","Names of the mesh adapter.");
   options.addConfigOption< std::vector<std::string> >("MeshAdapterMethod","Self-reg keys of the mesh adapter.");
   options.addConfigOption< std::vector<std::string> >("ErrorEstimatorNames","Names of the error estimator.");
   options.addConfigOption< std::vector<std::string> >("MeshCreator","Self-reg keys of the mesh creator.");
   options.addConfigOption< std::vector<std::string> >("ErrorEstimatorMethod","Self-reg keys of the error estimator.");
   options.addConfigOption< std::vector<std::string> >("OutputFormat","Self-reg keys of the output format.");
   options.addConfigOption< std::vector<std::string> >("MeshCreatorNames","Names of the mesh creator.");
   options.addConfigOption< std::vector<std::string> >("SpaceMethodNames","Names of the space method.");
   options.addConfigOption< std::vector<std::string> >("SpaceMethod","Self-reg keys of the space method.");
   options.addConfigOption< std::vector<std::string> >("OutputFormatNames","Names of the output format.");
   options.addConfigOption< CFuint >("InitialIter","Initial Iteration Number");
   options.addConfigOption< CFreal >("InitialTime","Initial Physical Time of the SubSystem");
}

//////////////////////////////////////////////////////////////////////////////

OnlyMeshSubSystem::OnlyMeshSubSystem(const std::string& name)
  : SubSystem(name),
    m_duration()
{
  addConfigOptionsTo(this);
  registActionListeners();

  // MeshCreator-related configuration options
  setParameter("MeshCreator",&m_meshCreator.mKeys);
  setParameter("MeshCreatorNames",&m_meshCreator.mNames);

  // MeshAdapter-related configuration options
  setParameter("MeshAdapterMethod",&m_meshAdapterMethod.mKeys);
  setParameter("MeshAdapterNames",&m_meshAdapterMethod.mNames);

  // ErrorEstimator-related configuration options
  setParameter("ErrorEstimatorMethod",&m_errorEstimatorMethod.mKeys);
  setParameter("ErrorEstimatorNames",&m_errorEstimatorMethod.mNames);

  // SpaceMethod-related configuration options
  setParameter("SpaceMethod",&m_spaceMethod.mKeys);
  setParameter("SpaceMethodNames",&m_spaceMethod.mNames);

  // OutputFormatter-related configuration options
  setParameter("OutputFormat",&m_outputFormat.mKeys);
  setParameter("OutputFormatNames",&m_outputFormat.mNames);

  m_initialTime = 0.;
  setParameter("InitialTime",&m_initialTime);

  m_initialIter = 0;
  setParameter("InitialIter",&m_initialIter);
}

//////////////////////////////////////////////////////////////////////////////

OnlyMeshSubSystem::~OnlyMeshSubSystem()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SubSystem::configure(args);

  // set the physical model
  configurePhysicalModel(args);

  // builds MeshCreator
  configureMultiMethod<MeshCreator>(args, m_meshCreator);

  // builds MeshAdapterMethod
  configureMultiMethod<MeshAdapterMethod>(args, m_meshAdapterMethod);

  // builds SpaceMethod
  configureMultiMethod<SpaceMethod>(args, m_spaceMethod);

  // builds ErrorEstimatorMethod
  configureMultiMethod<ErrorEstimatorMethod>(args, m_errorEstimatorMethod);

  // builds OutputFormatter
  configureMultiMethod<OutputFormatter>(args, m_outputFormat);

  // set collaborators
  setCollaborators();

  setParentNamespaceInMethodSockets(); // after this all the namespaces are set
  setEnableNamespaces(true);
  setNonRootMethods();
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::setNonRootMethods()
{
  CFAUTOTRACE;

  std::vector< Common::SafePtr<Method> > all_methods =
    MethodRegistry::getInstance().getAllMethods();

  std::for_each (all_methods.begin(),
                 all_methods.end(),
                 Common::safeptr_mem_fun(&Method::setNonRootMethods));
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::allocateAllSockets()
{
  CFAUTOTRACE;

  SubSystem::allocateAllSockets();

  // complement here the plugging of sockets and linking methods
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::deallocateAllSockets()
{
  CFAUTOTRACE;

  // complement here the unplugging of sockets and unlinking methods

  SubSystem::deallocateAllSockets();
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::buildPhysicalModel()
{
  CFAUTOTRACE;

  // first setup the physical models
  CFLog(NOTICE,"-------------------------------------------------------------\n");
  CFLog(NOTICE,"Setting up all PhysicalModel's\n\n");
  CFLog(NOTICE,"-------------------------------------------------------------\n");

  // loop on all namespaces and setup the models in each of them
  typedef std::vector<Common::SafePtr<Namespace> > NspVec;
  NspVec lst = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getAllNamespaces();
  cf_assert(!lst.empty()); // there should be some models
  
  NspVec::iterator nsp = lst.begin();
  for(; nsp != lst.end(); ++nsp)
  {
    NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).pushNamespace((*nsp)->getName());
    PhysicalModelStack::getActive()->getImplementor()->setup();
    NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).popNamespace();
  }
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::buildMeshData()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::buildMeshData() => START\n");
  
  // allocate all the mesh data
  vector <Common::SafePtr<MeshData> > meshDataVector = 
    MeshDataStack::getInstance().getAllEntries();

  for_each (meshDataVector.begin(), meshDataVector.end(), 
	    safeptr_mem_fun(&MeshData::reallocate));
    
  CFLog(NOTICE,"-------------------------------------------------------------\n");
  CFLog(NOTICE,"Setting up all MeshCreator's\n");
  CFLog(NOTICE,"-------------------------------------------------------------\n");
  
  if (getNbActiveRanks() < PE::GetPE().GetProcessorCount("Default")) {
    CFLog(INFO, "OnlyMeshSubSystem::buildMeshData() => [" << getNbActiveRanks() << "] active ranks\n"); 
  }
  
  if (PE::GetPE().GetRank("Default") < getNbActiveRanks()) {
    // finally setup the mesh creators
    Common::for_each_if (m_meshCreator.begin(), m_meshCreator.end(), 
			 mem_fun(&MeshCreator::setMethod),
			 mem_fun(&MeshCreator::isNonRootMethod), false);
    
    CFLog(NOTICE,"-------------------------------------------------------------\n");
    CFLog(NOTICE,"Building MeshData's\n");
    CFLog(NOTICE,"-------------------------------------------------------------\n");
    
    Common::for_each_if (m_meshCreator.begin(), m_meshCreator.end(), 
			 mem_fun(&MeshCreator::generateMeshData),
			 mem_fun(&MeshCreator::isNonRootMethod), false);
    
    // Process the CFMeshData to do:
    //  - renumbering
    //  - conversion FVM <-> FEM
    Common::for_each_if (m_meshCreator.begin(), m_meshCreator.end(), 
			 mem_fun(&MeshCreator::processMeshData),
			 mem_fun(&MeshCreator::isNonRootMethod), false);
    
    // Use the CFMeshData to build the mesh
   Common::for_each_if (m_meshCreator.begin(), m_meshCreator.end(), 
			mem_fun(&MeshCreator::buildMeshData),
			mem_fun(&MeshCreator::isNonRootMethod), false);
   
   // set up MeshData that are global across partitions 
   setGlobalMeshData();
   
   m_meshCreator.apply
     (root_mem_fun<void,MeshCreator>(&MeshCreator::unsetMethod));
  }
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::buildMeshData() => END\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::setup()
{
  CFAUTOTRACE;
 
  CFLog(VERBOSE, "OnlyMeshSubSystem::setup() start\n");
  
  cf_assert(isConfigured());
  
  if (PE::GetPE().GetRank("Default") < getNbActiveRanks()) {
    typedef std::vector<Common::SafePtr<Namespace> > NspVec;
    NspVec lst = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getAllNamespaces();
    cf_assert(!lst.empty());
    
    NspVec::iterator nsp = lst.begin();
    
    //setCommands() needs Trs's => must be exactly here
    setCommands();
    
    CFLog(NOTICE,"-------------------------------------------------------------\n");
    CFLogInfo("Setting up MeshAdapterMethod's\n");
    m_meshAdapterMethod.apply
      (mem_fun<void,MeshAdapterMethod>(&MeshAdapterMethod::setMethod));
    
    CFLog(NOTICE,"-------------------------------------------------------------\n");
    CFLogInfo("Setting up SpaceMethod's\n");
    m_spaceMethod.apply
      (mem_fun<void,SpaceMethod>(&SpaceMethod::setMethod));
    
    CFLog(NOTICE,"-------------------------------------------------------------\n");
    CFLogInfo("Setting up ErrorEstimatorMethod's\n");
    m_errorEstimatorMethod.apply
      (mem_fun<void,ErrorEstimatorMethod>(&ErrorEstimatorMethod::setMethod));
    
    CFLog(NOTICE,"-------------------------------------------------------------\n");
    CFLogInfo("Setting up OutputFormatter's\n");
    m_outputFormat.apply
      (root_mem_fun<void,OutputFormatter>(&OutputFormatter::setMethod));
    
    vector <Common::SafePtr<SubSystemStatus> > subSysStatusVec =
      SubSystemStatusStack::getInstance().getAllEntries();
    
    for(CFuint i = 0; i<subSysStatusVec.size();i++){
      // reset to 0 (or Initial Value) the number of iterations and the residual
      subSysStatusVec[i]->setNbIter(m_initialIter);
      subSysStatusVec[i]->setCurrentTimeDim(m_initialTime);
      // reset to 0 the number of iterations and the residual
      subSysStatusVec[i]->resetResidual();
    }
    
    // initialization has to be done after the creation of the normals
    // => after  m_spaceMethod->setMethod(), because some init
    // commands (init on the wall) need normals (init on the wall)
    CFLog(NOTICE,"\n-------------------------------------------------------------\n");
    CFLogInfo("Initializing solution\n");
    
    m_spaceMethod.apply(mem_fun<void,SpaceMethod>(&SpaceMethod::initializeSolution));
    
    CFLogInfo("Writing initial solution ... \n");
    writeSolution(true);
  }
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::setup() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::run()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "OnlyMeshSubSystem::run()\n");
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::unsetup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::unsetup() start\n");
  
  if (PE::GetPE().GetRank("Default") < getNbActiveRanks()) {
    Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
    
    vector <Common::SafePtr<SubSystemStatus> > subSysStatusVec =
      SubSystemStatusStack::getInstance().getAllEntries();
    
    // unset all the methods
    m_outputFormat.apply
      (root_mem_fun<void,OutputFormatter>(&OutputFormatter::unsetMethod));
    
    m_errorEstimatorMethod.apply
      (mem_fun<void,ErrorEstimatorMethod>(&ErrorEstimatorMethod::unsetMethod));
    
    m_spaceMethod.apply
      (mem_fun<void,SpaceMethod>(&SpaceMethod::unsetMethod));
    
    m_meshAdapterMethod.apply
      (mem_fun<void,MeshAdapterMethod>(&MeshAdapterMethod::unsetMethod));
  }
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::unsetup() end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

vector<Method*> OnlyMeshSubSystem::getMethodList()
{
  vector<Method*>  mList;

  for (CFuint i = 0; i < m_meshCreator.size(); ++i) {
    mList.push_back(m_meshCreator[i]);
  }

  for (CFuint i = 0; i < m_meshAdapterMethod.size(); ++i) {
    mList.push_back(m_meshAdapterMethod[i]);
  }

  for (CFuint i = 0; i < m_errorEstimatorMethod.size(); ++i) {
    mList.push_back(m_errorEstimatorMethod[i]);
  }

  for (CFuint i = 0; i < m_spaceMethod.size(); ++i) {
    mList.push_back(m_spaceMethod[i]);
  }

  for (CFuint i = 0; i < m_outputFormat.size(); ++i) {
    mList.push_back(m_outputFormat[i]);
  }

  return mList;
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::registActionListeners()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "OnlyMeshSubSystem::registActionListeners() start\n");
 
  // always call the parent class
  SubSystem::registActionListeners();

  // add here other ActionListeners
  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
  
  const std::string ssname = SubSystemStatusStack::getCurrentName();   
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MAESTRO_MODIFYRESTART"),
			     this,&OnlyMeshSubSystem::modifyRestartAction);
  event_handler->addListener(event_handler->key(ssname, "CF_ON_MESHADAPTER_AFTERGLOBALREMESHING"),
			     this,&OnlyMeshSubSystem::afterRemeshingAction);
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::registActionListeners() end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t OnlyMeshSubSystem::afterRemeshingAction(Common::Signal::arg_t eAfterRemesh)
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "OnlyMeshSubSystem::afterRemeshingAction() start\n");
 
  SimulationStatus::getInstance().setRestart(true);
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::afterRemeshingAction() end\n");
  
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t OnlyMeshSubSystem::modifyRestartAction(Common::Signal::arg_t eModifyRestart)
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "OnlyMeshSubSystem::modifyRestartAction() start\n");
  
  using namespace COOLFluiD::Environment;
  for(CFuint iMC=0; iMC < m_meshCreator.size(); iMC++)
  {
    if(!m_meshCreator[iMC]->isNonRootMethod())
    {
      std::string subsystemName = getName();
      std::string nspName = m_meshCreator[iMC]->getNamespace();
      std::string filename = SimulationStatus::getInstance().getLastOutputFile(subsystemName, nspName).string();
      CFLog(NOTICE,"Restarting from file : " << filename << "\n");

      //     CF_DEBUG_OBJ(Environment::DirPaths::getInstance().getWorkingDir().string());
      //     CF_DEBUG_OBJ(Environment::DirPaths::getInstance().getResultsDir().string());

      DirPaths::getInstance().setWorkingDir(
        DirPaths::getInstance().getResultsDir().string()
      );
      m_meshCreator[iMC]->modifyFileNameForRestart(filename);
    }
  }

  m_initialTime = SimulationStatus::getInstance().getSimulationTime(getName());
  vector< SafePtr<SubSystemStatus> > subSystemStatusVec =
    SubSystemStatusStack::getInstance().getAllEntries();

  for(CFuint i=0; i<subSystemStatusVec.size(); ++i)
  {
    subSystemStatusVec[i]->setCurrentTimeDim(m_initialTime);
  }
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::modifyRestartAction() end\n");
  
  return Common::Signal::return_t ();
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::setCollaborators()
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"Setting collaborators for all Methods \n");

  CFLog(DEBUG_MIN,"Setting OutputFormatter collaborators for MeshAdapterMethods \n");
  SubSystem::setCollaborators<MeshAdapterMethod, OutputFormatter>(m_meshAdapterMethod, m_outputFormat);

  CFLog(DEBUG_MIN,"Setting MeshCreator collaborators for MeshAdapterMethods \n");
  SubSystem::setCollaborators<MeshAdapterMethod, MeshCreator>(m_meshAdapterMethod, m_meshCreator);

  CFLog(DEBUG_MIN,"Setting SpaceMethod collaborators for MeshAdapterMethods \n");
  SubSystem::setCollaborators<MeshAdapterMethod, SpaceMethod>(m_meshAdapterMethod, m_spaceMethod);
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::writeSolution(const bool force_write )
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "OnlyMeshSubSystem::writeSolution() start\n");
  
  for (unsigned int i = 0; i < m_outputFormat.size(); ++i)
  {
    if(!m_outputFormat[i]->isNonRootMethod())
    {
      if ( m_outputFormat[i]->isSaveNow( force_write ) )
      {
        Stopwatch<WallTime> stopTimer;
        stopTimer.start();
        m_outputFormat[i]->open ();
        m_outputFormat[i]->write();
        m_outputFormat[i]->close();
        stopTimer.stop();
        CFLog(INFO, "Writing took " << stopTimer << "s\n");
      }
    }
  }
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::writeSolution() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::configurePhysicalModel ( Config::ConfigArgs& args )
{
  CFLog(VERBOSE, "OnlyMeshSubSystem::configurePhysicalModel() start\n");
 
  typedef std::vector<Common::SafePtr<Namespace> > NspVec;
  NspVec lst = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getAllNamespaces();
  
  cf_assert(!lst.empty());

  NspVec::iterator nsp = lst.begin();

  for(; nsp != lst.end(); ++nsp)
  {
    const std::string physicalModelName = (*nsp)->getPhysicalModelName();
    const std::string physicalModelType = (*nsp)->getPhysicalModelType();

    // create the new physical model implementor
    Common::SafePtr<PhysicalModelImpl::PROVIDER> physicalMdlProv =
      FACTORY_GET_PROVIDER(getFactoryRegistry(), PhysicalModelImpl, physicalModelType);
    cf_assert(physicalMdlProv.isNotNull());
    
    Common::SelfRegistPtr<PhysicalModelImpl> physicalModelImpl =
      physicalMdlProv->create(physicalModelName);
    cf_assert(physicalModelImpl.isNotNull());
    
    physicalModelImpl->setFactoryRegistry(getFactoryRegistry());
    
    // configure the physical model implementor
    configureNested ( physicalModelImpl.getPtr(), args );

    // set the PhysicalModelImpl in the PhysicalModel
    // this phase has to be done after the configuration because
    // some data (example: nbEquations for models with multiple species)
    // could have required a configuration phase first
    PhysicalModelStack::getInstance().getEntry(physicalModelName)->setPhysicalModelImpl(physicalModelImpl);

    // default initialization of the equation subsystem descriptor
    const CFuint nbEqs = PhysicalModelStack::getInstance().getEntry(physicalModelName)->getNbEq();
    PhysicalModelStack::getInstance().getEntry(physicalModelName)->setEquationSubSysDescriptor(0,nbEqs,0);
  }
  
  CFLog(VERBOSE, "OnlyMeshSubSystem::configurePhysicalModel() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::writeConvergenceOnScreen()
{
}

//////////////////////////////////////////////////////////////////////////////

void OnlyMeshSubSystem::setGlobalMeshData()
{
  vector <Common::SafePtr<MeshData> > meshDataVector = 
    MeshDataStack::getInstance().getAllEntries();
  
  const int rank = PE::GetPE().GetRank("Default");
  for(CFuint iMeshData = 0; iMeshData < meshDataVector.size(); iMeshData++) {
    SafePtr<MeshData> currMeshData = meshDataVector[iMeshData];
    bool isNonRoot = false;
    for(CFuint iMC = 0; iMC < m_meshCreator.size(); iMC++) {
      if(m_meshCreator[iMC]->getNamespace() == currMeshData->getPrimaryNamespace())
        isNonRoot = m_meshCreator[iMC]->isNonRootMethod();
    }
    
    if (!isNonRoot) {
      // AL: does this work if the same mesh is shared by different namespaces??
      if (PE::GetPE().isRankInGroup(rank, currMeshData->getPrimaryNamespace())) {
	CFLogNotice("Building TRS info for MeshData in Namespace " << currMeshData->getPrimaryNamespace() << "\n");
	// === conversion fix ===
	vector<string>& TotalTRSNames = currMeshData->getTotalTRSNames ();
	vector<vector<CFuint> >& TotalTRSInfo = currMeshData->getTotalTRSInfo ();
	vector< SafePtr<TopologicalRegionSet> > trsList = currMeshData->getTrsList();
	
	if (TotalTRSInfo.empty ()) {
	  cf_assert(TotalTRSNames.empty());
	  
	  // count in advance the number of writable TRS
	  CFuint sizeTRSToWrite = 0;
	  for (CFuint i = 0; i < trsList.size(); ++i) {
	    if (trsList[i]->hasTag("writable")) {
	      sizeTRSToWrite++;
	    }
	  }
	  cf_assert(sizeTRSToWrite <= trsList.size());
	  
	  TotalTRSNames.resize(sizeTRSToWrite);
	  TotalTRSInfo.resize(sizeTRSToWrite);
	  
	  CFuint counter = 0;
	  vector< SafePtr<TopologicalRegionSet> >::iterator it;
	  for (it = trsList.begin(); it != trsList.end(); ++it) {
	    if ((*it)->hasTag("writable")) {
	      TotalTRSNames[counter] = (*it)->getName();
	      const CFuint nbTRs = (*it)->getNbTRs();
	      TotalTRSInfo[counter].resize(nbTRs);
	      for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
		// we make sure that the number of boundary faces is always updated
		TotalTRSInfo[counter][iTR] = (*(*it))[iTR]->getLocalNbGeoEnts();
		CFLog(VERBOSE, "TRS : " << (*it)->getName()
		      << ", TR : "<< iTR
		      << ", nbGeos : "
		      << TotalTRSInfo[counter][iTR] << "\n");
	      }
	      counter++;
	    }
	  }
	  
	  cf_assert(counter == sizeTRSToWrite);
	}
	
	if (PE::GetPE().IsParallel() && 
	    (currMeshData->getNbNodes() > 0 && currMeshData->getNbStates() > 0)) {
	  const std::string parStateVecName = currMeshData->getPrimaryNamespace() + "_states";
	  const std::string parNodeVecName = currMeshData->getPrimaryNamespace() + "_nodes";
	  
	  DataHandle<State*, GLOBAL> states =
	    currMeshData->getDataStorage()->getGlobalData<State*>(parStateVecName);
	  
	  DataHandle<Node*, GLOBAL> nodes =
	    currMeshData->getDataStorage()->getGlobalData<Node*>(parNodeVecName);
	  
	  // AL : this should be activated only in very special occasions
	  // #ifndef NDEBUG
	  //  states.DumpContents ();
	  //  nodes.DumpContents ();
	  //  #endif
	  
	  states.buildMap(CFEnv::getInstance().getVars()->SyncAlgo);
	  nodes.buildMap(CFEnv::getInstance().getVars()->SyncAlgo);
	  
	  // #ifndef NDEBUG
	  //  states.DumpInfo ();
	  //  nodes.DumpInfo ();
	  // #endif
	}
      }
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
