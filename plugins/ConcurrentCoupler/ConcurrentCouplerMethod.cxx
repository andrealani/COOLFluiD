#include "Common/FilesystemException.hh"
#include "Common/BadValueException.hh"
#include "Common/PEFunctions.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/DirPaths.hh"

#include "Framework/CouplerMethod.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "ConcurrentCoupler/ConcurrentCoupler.hh"
#include "ConcurrentCoupler/ConcurrentCouplerMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ConcurrentCouplerMethod,
			    CouplerMethod,
			    ConcurrentCouplerModule,
			    1>
concurrentCouplerMethodProvider("ConcurrentCoupler");

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");
  
  options.addConfigOption< vector<string> >("InterfacesReadNames","Names for the configuration of the commands.");
  options.addConfigOption< vector<string> >("InterfacesReadComs","dataTransferRead commands.");
 
  options.addConfigOption< vector<string> >("InterfacesWriteNames","Names for the configuration of the commands.");
  options.addConfigOption< vector<string> >("InterfacesWriteComs","dataTransferWrite commands.");
  
  options.addConfigOption< vector<string> >("CoupledSubSystems","Names of the subsystems to be coupled to.");
  options.addConfigOption< vector<string> >("CoupledNameSpaces","Names of the namespaces of the other subsystems.");
  
  options.addConfigOption< vector<CFuint> >("TransferRates","Transfer data every X iterations for each namespace.");
}
      
//////////////////////////////////////////////////////////////////////////////

ConcurrentCouplerMethod::ConcurrentCouplerMethod(const string& name) :
  CouplerMethod(name),
  m_setup(),
  m_interfacesRead(CFNULL),
  m_interfacesWrite(CFNULL),
  m_fileRW(CFNULL, CFNULL),
  m_ioRoot(false)
{
  addConfigOptionsTo(this);
  
  m_data.reset(new ConcurrentCouplerData(this));
  
  m_setupStr = "Null";
  setParameter("SetupCom",&m_setupStr);
  
  m_interfacesReadStr = vector<string>();
  setParameter("InterfacesReadComs",&m_interfacesReadStr);
  m_interfacesReadNameStr = vector<string>();
  setParameter("InterfacesReadNames",&m_interfacesReadNameStr);
  
  m_interfacesWriteStr = vector<string>();
  setParameter("InterfacesWriteComs",&m_interfacesWriteStr);
  m_interfacesWriteNameStr = vector<string>();
  setParameter("InterfacesWriteNames",&m_interfacesWriteNameStr);
  
  m_coupledSubSystemsStr = vector<string>();
  setParameter("CoupledSubSystems",&m_coupledSubSystemsStr);
  m_coupledNamespacesStr = vector<string>();
  setParameter("CoupledNameSpaces",&m_coupledNamespacesStr);
  
  m_transferRates = vector<CFuint>();
  setParameter("TransferRates",&m_transferRates);
}
      
//////////////////////////////////////////////////////////////////////////////

ConcurrentCouplerMethod::~ConcurrentCouplerMethod()
{
  clearComds();
  
  delete m_fileRW.first;
  delete m_fileRW.second;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> ConcurrentCouplerMethod::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::configure ( Config::ConfigArgs& args )
{
  CouplerMethod::configure(args);
    
  cf_assert (m_coupledNamespacesStr.size() > 0);
  
  if (m_coupledSubSystemsStr.size() == 0) {
    const string ssname = SubSystemStatusStack::getCurrentName();
    m_coupledSubSystemsStr.resize(m_coupledNamespacesStr.size(), ssname);
  }
  
  cf_assert (m_coupledNamespacesStr.size() == m_coupledSubSystemsStr.size());
  
  // Check size of configuration parameters
  // cf_assert(m_interfacesReadStr.size()    == m_interfacesReadNameStr.size());
  cf_assert(m_interfacesWriteStr.size()   == m_interfacesWriteNameStr.size());
    
  configureNested ( m_data.getPtr(), args );
  
  // configureInterfaces();
  
  // add configures to the ConcurrentCouplerCom's
  clearComds();

  // Make the name accessible to all commands by putting them in the data
  // m_data->setCoupledSubSystemsNames(m_coupledSubSystemsStr, m_coupledNamespacesStr);
  // m_data->resizeCoupledInterfaces();
  
  // // Resize the vectors of commands
  if (m_interfacesReadStr.size() > 0) {
    m_interfacesRead.resize(m_interfacesReadStr.size());
  }
  
  if (m_interfacesWriteStr.size() > 0) {
    m_interfacesWrite.resize(m_interfacesWriteStr.size());
  }
  
  if (m_transferRates.size() != m_coupledNamespacesStr.size()) {
    if (m_transferRates.size() == 0) {
      CFLog(WARN, "ConcurrentCouplerMethod::configure() => m_transferRates=1 for all namespaces\n");
      m_transferRates.resize(m_coupledNamespacesStr.size(),1);
    }
    else {
      throw BadValueException
	(FromHere(), "ConcurrentCouplerMethod::configure() => m_transferRates.size() != m_coupledNamespacesStr.size()\n");
    }
  }
  
  for (CFuint i = 0; i < m_transferRates.size(); ++i) {
    cf_assert(m_transferRates[i] > 0);
  }
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ConcurrentCouplerMethod::postConfigure( Config::ConfigArgs& args )
{
  configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
    (args, m_setup,m_setupStr,m_setupStr,m_data);
  
  CFLog(VERBOSE, "ConcurrentCoupler [" << getName() << "] Namespace [" << getNamespace()
	<< "] Data namespace [" << getMethodData()->getNamespace() << "]\n");
  
  // Configure all commands
  for (CFuint i = 0; i < m_interfacesReadStr.size(); ++i) {
    CFLog(INFO, "Interface Reader type = " << m_interfacesReadStr[i] << "\n");
    CFLog(INFO, "Configuring Interface Reader command ["  << m_interfacesReadStr[i]   << "]\n");
    configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
      (args,m_interfacesRead[i],m_interfacesReadStr[i], m_interfacesReadNameStr[i], m_data);
    cf_assert(m_interfacesRead[i].isNotNull());
  }
  
  for (CFuint i = 0; i < m_interfacesWriteStr.size(); ++i) {
    CFLog(INFO, "Interface Writer type = " << m_interfacesWriteStr[i] << "\n");
    CFLog(INFO, "Configuring Interface Writer command ["  << m_interfacesWriteStr[i]   << "]\n");
    configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
      (args,m_interfacesWrite[i],m_interfacesWriteStr[i], m_interfacesWriteNameStr[i], m_data);
    cf_assert(m_interfacesWrite[i].isNotNull());
  }
  
  // // Check that all the commands apply to the same group of TRSs
  // for(CFuint i = 0; i < m_interfacesRead.size(); ++i) {
  //   cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_matchMeshesRead[i]->getCommandGroupName());
  //   cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_matchMeshesWrite[i]->getCommandGroupName());
  //   cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_preProcessRead[i]->getCommandGroupName());
  //   cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_preProcessWrite[i]->getCommandGroupName());
  //   cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_postProcess[i]->getCommandGroupName());
  //   cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_matchMeshesRead[i]->getCommandGroupName());
  //   cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_matchMeshesWrite[i]->getCommandGroupName());
  //   cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_preProcessRead[i]->getCommandGroupName());
  //   cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_preProcessWrite[i]->getCommandGroupName());
  //   cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_postProcess[i]->getCommandGroupName());
  // }
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::configureInterfaces()
{
  // All Command Groups are assumed to be interfaces
  // m_data->setInterfaces(getCommandGroups());
  
  // // prepare the InterfacesMap vector
  // SafePtr<CFMap<string, CFuint> > coupledInterfacesMap = m_data->getCoupledInterfacesMapPtr();
  // coupledInterfacesMap->reserve(m_groups.size());
  
  // // we need to have at least one group
  // cf_assert(getNbGroups() > 0);

  // // fill the map
  // for(CFuint i = 0; i < m_groups.size(); ++i) {
  //   coupledInterfacesMap->insert(m_groups[i]->getName(),i);
  // }
  
  // coupledInterfacesMap->sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::setInfoToOtherSubSystem()
{
  CFAUTOTRACE;
  
  const string nameSpace = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(nameSpace);
  
  // /// @todo this will be substituted by a communication stream
  // const bool isParallel = Common::PE::GetPE().IsParallel();
  // if((!isParallel) || ((isParallel) && (Common::PE::GetPE().GetRank(nameSpace) == 0)))
  // {
  //   Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
  //   Common::SafePtr<PhysicalModel> physicalModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  //   const string currentSubSystem = subsystemStatus->getSubSystemName();

  //   // For each interface, write the number and names of TRS
  //   for(CFuint i = 0; i < m_groups.size(); ++i) {

  //     const string interfaceName = m_groups[i]->getName();

  //     const vector<string>& trsNames = m_groups[i]->getTrsNames();
  //     const CFuint nbTRS = trsNames.size();

  //     boost::filesystem::path dataHandleName ("COUPLING_" + interfaceName + "_" + nameSpace + "_" + currentSubSystem);
  //     boost::filesystem::path nameOutputFile =
  //       Environment::DirPaths::getInstance().getResultsDir() / dataHandleName;

  //     SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  //     ofstream& fout = fhandle->open(nameOutputFile);
      
  //     const vector<string> coordType = m_data->getThisCoordType(interfaceName);
  //     //Get variable transformer
  //     Common::SafePtr<PreVariableTransformer> varTransfo = m_data->getPreVariableTransformer(interfaceName);
      
  //     fout << "!COORDTYPE " << coordType.size() << "\n";
  //     for (CFuint i = 0; i < coordType.size(); ++i) {
  //       fout << coordType[i] << "\n";
  //     }
  //     fout << "!DATASIZE " << varTransfo->getTransformedSize(physicalModel->getNbEq()) << "\n";
      
  //     fout << "!NBTRS " << nbTRS << "\n";
  //     for (CFuint i = 0; i < nbTRS; ++i) {
  //       fout << trsNames[i] << "\n";
  //     }
      
  //     fhandle->close();
  //   }
  // }
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::getInfoFromOtherSubSystem()
{
  CFAUTOTRACE;
  
  // const string nsp = getMethodData()->getNamespace();  
  // runSerial<void, ConcurrentCoupler, &ConcurrentCouplerMethod::readInfoFromOtherSubSystem>(this, nsp);
}
      
//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::readInfoFromOtherSubSystem()
{

  // Get the information about the other subsystem TRS's coupled
  // for(CFuint i = 0; i < m_groups.size(); ++i) {

  //   const string interfaceName = m_groups[i]->getName();

  //   const string otherNamespace = m_data->getCoupledNameSpaceName(interfaceName);
  //   const string otherSubSystem = m_data->getCoupledSubSystemName(interfaceName);
  //   boost::filesystem::path fileName ("COUPLING_" + interfaceName + "_" + otherNamespace + "_" + otherSubSystem);

  //   //Reading the file
  //   boost::filesystem::path fname =
  //     Environment::DirPaths::getInstance().getResultsDir() / fileName;
    
  //   Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = 
  //     Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  //   ifstream& fin = fhandle->open(fname);
    
  //   CFuint lineNb = 0;
  //   string line;
  //   vector<string> words;

  //   // read type of Data
  //   getWordsFromLine(fin,line,lineNb,words);
  //   cf_assert(words.size() == 2);
  //   cf_assert(words[0] == "!COORDTYPE");
  //   const CFuint nbCoordType = Common::StringOps::from_str<CFint>(words[1]);
  //   vector<string> coordType(nbCoordType);
  //   for (CFuint i=0; i< nbCoordType; ++i)
  //   {
  //     getWordsFromLine(fin,line,lineNb,words);
  //     cf_assert(words.size() == 1);
  //     coordType[i] = words[0];
  //     cf_assert((coordType[i] == "Nodes")||(coordType[i] == "States") ||(coordType[i] == "Ghost") ||(coordType[i] == "Gauss"));
  //   }

  //   // read size of Data
  //   getWordsFromLine(fin,line,lineNb,words);
  //   cf_assert(words.size() == 2);
  //   cf_assert(words[0] == "!DATASIZE");
  //   const CFuint dataSize = Common::StringOps::from_str<CFint>(words[1]);

  //   // read nb trs's
  //   getWordsFromLine(fin,line,lineNb,words);
  //   cf_assert(words.size() == 2);
  //   cf_assert(words[0] == "!NBTRS");
  //   const CFuint nbTRS = Common::StringOps::from_str<CFint>(words[1]);
  //   vector<string> trsNames(nbTRS);

  //   for (CFuint i=0; i< nbTRS; ++i)
  //   {
  //     getWordsFromLine(fin,line,lineNb,words);
  //     cf_assert(words.size() == 1);
  //     trsNames[i] = words[0];
  //   }

  //   fhandle->close();

  //   m_data->setCoupledSubSystemTRSNames(interfaceName, trsNames);
  //   m_data->setCoupledSubSystemCoordTypes(interfaceName, coordType);
  //   m_data->setTransferedSize(interfaceName, dataSize);
  // }
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::setMethodImpl()
{
  CFAUTOTRACE;

  CouplerMethod::setMethodImpl();
  
  //first check that the method and its commands have been correctly configured
  cf_assert(isConfigured());
  cf_assert(m_setup.isNotNull());
  
  m_setup->execute();
  
  setupCommandsAndStrategies(); 
  
  const string fileName = getStatusFilename();
  if (PE::GetPE().GetRank(getNamespace()) == 0) {
    ofstream fout(fileName.c_str());
    // initialize the file with all 0's
    resetStatusFile(&fout, READ);
    resetStatusFile(&fout, WRITE);
    fout.close();
  }
  
  Group& cgroup = PE::GetPE().getGroup(getNamespace());
  MPI_Barrier(cgroup.comm);
  
  // loop over all coupled namespaces to find a root I/O for each one of them 
  const int rank  = PE::GetPE().GetRank("Default");     // rank in default group
  const CFuint nspSize = m_coupledNamespacesStr.size();
  for (CFuint i = 0; i < nspSize; ++i) {
    if (PE::GetPE().isRankInGroup(rank, m_coupledNamespacesStr[i])) { 
      // select process with group rank=0 as root I/O
      if (PE::GetPE().GetRank(m_coupledNamespacesStr[i]) == 0) {
	m_ioRoot = true;
	CFLog(INFO, "ConcurrentCouplerMethod::setMethodImpl() => namespace [" << 
	      m_coupledNamespacesStr[i] << "] has root I/O rank [" << rank << "]\n");
	
	m_fileRW.second = new ofstream(fileName.c_str(), ios_base::in | ios_base::out);
	m_fileRW.second->close();
	m_fileRW.first  = new ifstream(fileName.c_str(), ios_base::in);	
	m_fileRW.first->close();
      } 
      break;
    }
  } 
}
      
//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::preProcessReadImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());
  
  // if(Common::PE::GetPE().IsParallel()){
  //   Common::PE::GetPE().setBarrier(getMethodData()->getNamespace());
  // }
  
  // // preprocess interfaces
  // for(CFuint i = 0; i < m_preProcessRead.size(); ++i) {
  //   cf_assert(m_preProcessRead[i].isNotNull());
  //   m_preProcessRead[i]->execute();
  // }
  
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::preProcessWriteImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());

  // // preprocess interfaces
  // for(CFuint i = 0; i < m_preProcessWrite.size(); ++i) {
  //   cf_assert(m_preProcessWrite[i].isNotNull());
  //   m_preProcessWrite[i]->execute();
  // }

}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::meshMatchingReadImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());
  
  // const string nsp = getMethodData()->getNamespace();
  // if (Common::PE::GetPE().IsParallel())
  // {
  //   Common::PE::GetPE().setBarrier(nsp);
  // }

  // // match interface meshes
  // for(CFuint i = 0; i < m_matchMeshesRead.size(); ++i) {
  //   cf_assert(m_matchMeshesRead[i].isNotNull());
  //   m_matchMeshesRead[i]->execute();
  // }

}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::meshMatchingWriteImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());

  // // match interface meshes
  // for(CFuint i = 0; i < m_matchMeshesWrite.size(); ++i) {
  //   cf_assert(m_matchMeshesWrite[i].isNotNull());
  //   m_matchMeshesWrite[i]->execute();
  // }
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::dataTransferReadImpl()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "ConcurrentCouplerMethod::dataTransferReadImpl() => start\n");
  
  cf_assert(isSetup());
  cf_assert(isConfigured());
  
  //  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const bool doCoupling = isCouplingIter(m_fileRW, READ);
  if (doCoupling) {
    for(CFuint i = 0; i < m_interfacesRead.size(); ++i) {
      cf_assert(m_interfacesRead[i].isNotNull());
      // Execute and save file if needed...
      //    if((!(iter % m_transferRates[i])) || (iter ==0) ) {
      m_interfacesRead[i]->execute();
    }
    
    if (PE::GetPE().GetRank(getNamespace()) == 0) {
      const string fileName = getStatusFilename();
      ofstream fout(fileName.c_str());
      // initialize the file with all 0's
      resetStatusFile(&fout, READ);
    } 
    
    Group& cgroup = PE::GetPE().getGroup(getNamespace());
    MPI_Barrier(cgroup.comm);
  }
  
  CFLog(INFO, "ConcurrentCouplerMethod::dataTransferReadImpl() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::dataTransferWriteImpl()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "ConcurrentCouplerMethod::dataTransferWriteImpl() => start\n");
  
  cf_assert(isSetup());
  cf_assert(isConfigured());
  
  const bool doCoupling = isCouplingIter(m_fileRW, WRITE);
  if (doCoupling) {
    CFLog(INFO, "ConcurrentCouplerMethod::dataTransferWriteImpl() => m_interfacesWrite.size() = " 
	  << m_interfacesWrite.size() << "\n");
    
    for(CFuint i = 0; i < m_interfacesWrite.size(); ++i) {
      cf_assert(m_interfacesWrite[i].isNotNull());
      
      CFLog(INFO, "ConcurrentCouplerMethod::dataTransferWriteImpl() => before [" 
	    << m_interfacesWrite[i]->getName() << "]\n");
      m_interfacesWrite[i]->execute();
      CFLog(INFO, "ConcurrentCouplerMethod::dataTransferWriteImpl() => after [" 
	    << m_interfacesWrite[i]->getName() << "]\n");
      
      // barrier ensures that at this point all coupled ranks synchronize
      CFLog(INFO, "before barrier\n");
      Group& cgroup = PE::GetPE().getGroup(getNamespace());
      MPI_Barrier(cgroup.comm);
      CFLog(INFO, "after barrier\n");
    }
    
    if (PE::GetPE().GetRank(getNamespace()) == 0) {
      const string fileName = getStatusFilename();
      ofstream fout(fileName.c_str());
      // initialize the file with all 0's
      resetStatusFile(&fout, WRITE);
    }   
  }
  
  CFLog(INFO, "ConcurrentCouplerMethod::dataTransferWriteImpl() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::unsetMethodImpl()
{
  unsetupCommandsAndStrategies();

  CouplerMethod::unsetMethodImpl();
}
      
//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::clearComds()
{
  CFAUTOTRACE;
  
  if (m_interfacesRead.size() > 0) {
    vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_interfacesRead);
  }
  
  if (m_interfacesWrite.size() > 0) {
    vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_interfacesWrite);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<Common::SafePtr<Framework::NumericalStrategy> >
ConcurrentCouplerMethod::getStrategyList() const
{
  vector<Common::SafePtr<Framework::NumericalStrategy> > result;

  // vector < Common::SelfRegistPtr<PreVariableTransformer > > preTrans = m_data->getPreVariableTransformers();
  // for(CFuint i = 0; i < preTrans.size(); ++i) {
  //   result.push_back(preTrans[i].getPtr());
  // }
  
  // vector < Common::SelfRegistPtr<PostVariableTransformer > > postTrans = m_data->getPostVariableTransformers();
  // for(CFuint i = 0; i < postTrans.size(); ++i) {
  //   result.push_back(postTrans[i].getPtr());
  // }
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::getWordsFromLine(ifstream& fin,
					       string& line,
					       CFuint& lineNb,
					       vector<string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

bool ConcurrentCouplerMethod::isCouplingIter(pair<ifstream*, ofstream*>& file, 
					     TypeIO tio)
{
  cf_assert(m_transferRates.size() == m_coupledNamespacesStr.size());
  
  // There is a ASCII file named "status.<CurrentNamespace>", where <CurrentNamespace> 
  // is returned by ConcurrentCouplerMethod::getNamespace().
  // In this file a flag (0 or 1) is associated to each coupled namespace (the ones 
  // user-defined with the "CoupledNameSpaces" option) within <CurrentNamespace> for 
  // reading (ConcurrentCouplerMethod::dataTransferReadImpl()) and 
  // writing (ConcurrentCouplerMethod::dataTransferWriteImpl()) coupling phases 
  //
  // Examples: case with 3 coupled namespaces (3 flags for reading, 3 flags for writing)
  // 010010 => only the second namespace is ready (1), so no coupling will occur
  // 111000 => all 3 namespaces are ready for coupling (only reading phase, no writing)
  // 000111 => all 3 namespaces are ready for coupling (only writing phase, no reading)
  // 111111 => all 3 namespaces are ready for coupling (reading and writing phases)
  //
  // In each group corresponding to a coupled namespace, the root process will
  // 1) write only the flag value for its group in the corresponding position
  // 2) read at once all flags in the file
  // 3) communicate a flag telling whether to couple to all other ranks in its group
  
  // loop over all coupled namespaces to find the one containing the current rank
  const int rank  = PE::GetPE().GetRank("Default");     // rank in default group
  const CFuint nspSize = m_coupledNamespacesStr.size();
  bool ready = false;
  int nspID = -1;
  vector<char> flags(nspSize);
  const string fileName = getStatusFilename();	
  
  for (CFuint i = 0; i < nspSize; ++i) {
    if (PE::GetPE().isRankInGroup(rank, m_coupledNamespacesStr[i])) { 
      SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
	(SubSystemStatusStack::getCurrentName()).getNamespace(m_coupledNamespacesStr[i]);
      SafePtr<SubSystemStatus> subSysStatus = 
	SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
      cf_assert(subSysStatus.isNotNull());
      nspID = i;
      // store the number of iterations in each of the namespaces (to be coupled)
      // this is =0 f the namespace is locally not active
      cf_assert(m_transferRates[i] > 0);
      const CFuint nbIter = subSysStatus->getNbIter();
      if (nbIter%m_transferRates[i] == 0) {
	ready = true;
	if (m_ioRoot) {
	  const long offset = tio*nspSize*sizeof(bool);
	  file.second->open(fileName.c_str(), ios_base::in | ios_base::out);
	  file.second->seekp(i*sizeof(bool) + offset);
	  (*file.second) << ready;
	  file.second->close();
	}
	break;
      }
    }
  }
  
  cf_assert(nspID >= 0);
    
  int doCoupling = 1;
  if (m_ioRoot) {
    const long offset = tio*nspSize*sizeof(bool);
    file.first->open(fileName.c_str(), ios_base::in);
    file.first->seekg(offset);
    for (CFuint i = 0; i < nspSize; ++i) {
      (*file.first) >> flags[i];
      CFLog(VERBOSE, "ConcurrentCouplerMethod::isCouplingIter() on P[" << rank 
	    << "] at offset [" << offset << "] => flag [" << i << "]=" << flags[i] << "\n");
      // if at least one flag is not "true" then get out
      if (flags[i] == '0') {
	doCoupling = 0;
	break;
      }
    }
    file.first->close();
  }
  
  // Communicate to all the processors belonging to the same coupled namespace 
  Group& group = PE::GetPE().getGroup(m_coupledNamespacesStr[nspID]);
  MPIError::getInstance().check
    ("MPI_Bcast", "ConcurrentCouplerMethod::isCouplingIter()", 
     MPI_Bcast(&doCoupling, 1, MPI_INT, 0, group.comm));
  
  CFLog(VERBOSE, "ConcurrentCouplerMethod::isCouplingIter() on P[" << rank 
	<< "] is ready for coupling? " << doCoupling << "\n");
  
  return doCoupling;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler 

  } // namespace Numerics

} // namespace COOLFluiD

