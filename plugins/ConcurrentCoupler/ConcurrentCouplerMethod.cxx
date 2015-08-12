#include "Common/FilesystemException.hh"
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
  
  options.addConfigOption< vector<string> >("PreProcessReadComs","Preprocessing");
  options.addConfigOption< vector<string> >("PreProcessReadNames","Names for the configuration of the commands.");
  options.addConfigOption< vector<string> >("PreProcessWriteComs","Preprocessing");
  options.addConfigOption< vector<string> >("PreProcessWriteNames","Names for the configuration of the commands.");
  
  options.addConfigOption< vector<string> >("MeshMatchingReadComs","MeshMatching command.");
  options.addConfigOption< vector<string> >("MeshMatchingReadNames","Names for the configuration of the commands.");
  options.addConfigOption< vector<string> >("MeshMatchingWriteComs","MeshMatching command.");
  options.addConfigOption< vector<string> >("MeshMatchingWriteNames","Names for the configuration of the commands.");
  
  options.addConfigOption< vector<string> >("InterfacesReadNames","Names for the configuration of the commands.");
  options.addConfigOption< vector<string> >("InterfacesReadComs","dataTransferRead commands.");
  options.addConfigOption< vector<string> >("InterfacesWriteNames","Names for the configuration of the commands.");
  options.addConfigOption< vector<string> >("InterfacesWriteComs","dataTransferWrite commands.");
  
  options.addConfigOption< vector<string> >("PostProcessNames","Names for the configuration of the commands.");
  options.addConfigOption< vector<string> >("PostProcessComs","Unsetup commands.");
  
  options.addConfigOption< vector<string> >("InterfacesNames","Names of the interfaces.");
  options.addConfigOption< vector<string> >("CoupledSubSystems","Names of the subsystems to be coupled to.");
  options.addConfigOption< vector<string> >("CoupledNameSpaces","Names of the namespaces of the other subsystems.");
  
  options.addConfigOption< vector<CFuint> >("TransferRates","Transfer data every X iterations");
}
      
//////////////////////////////////////////////////////////////////////////////

ConcurrentCouplerMethod::ConcurrentCouplerMethod(const string& name) :
  CouplerMethod(name),
  m_preProcessRead(0),
  m_preProcessWrite(0),
  m_matchMeshesRead(0),
  m_matchMeshesWrite(0),
  m_interfacesRead(0),
  m_interfacesWrite(0),
  m_postProcess(0),
  m_fullConfigure(0)
{
  addConfigOptionsTo(this);
  m_data.reset(new ConcurrentCouplerData(this));
    
  m_setupStr = "Null";
  setParameter("SetupCom",&m_setupStr);
  
  m_preProcessReadStr = vector<string>();
  setParameter("PreProcessReadComs",&m_preProcessReadStr);
  m_preProcessReadNameStr = vector<string>();
  setParameter("PreProcessReadNames",&m_preProcessReadNameStr);
  
  m_preProcessWriteStr = vector<string>();
  setParameter("PreProcessWriteComs",&m_preProcessWriteStr);
  m_preProcessWriteNameStr = vector<string>();
  setParameter("PreProcessWriteNames",&m_preProcessWriteNameStr);
  
  m_matchMeshesReadStr = vector<string>();
  setParameter("MeshMatchingReadComs",&m_matchMeshesReadStr);
  m_matchMeshesReadNameStr = vector<string>();
  setParameter("MeshMatchingReadNames",&m_matchMeshesReadNameStr);
  
  m_matchMeshesWriteStr = vector<string>();
  setParameter("MeshMatchingWriteComs",&m_matchMeshesWriteStr);
  m_matchMeshesWriteNameStr = vector<string>();
  setParameter("MeshMatchingWriteNames",&m_matchMeshesWriteNameStr);
  
  m_postProcessStr = vector<string>();
  setParameter("PostProcessComs",&m_postProcessStr);
  m_postProcessNameStr = vector<string>();
  setParameter("PostProcessNames",&m_postProcessNameStr);
  
  m_interfacesReadStr = vector<string>();
  setParameter("InterfacesReadComs",&m_interfacesReadStr);
  m_interfacesReadNameStr = vector<string>();
  setParameter("InterfacesReadNames",&m_interfacesReadNameStr);
  
  m_interfacesWriteStr = vector<string>();
  setParameter("InterfacesWriteComs",&m_interfacesWriteStr);
  m_interfacesWriteNameStr = vector<string>();
  setParameter("InterfacesWriteNames",&m_interfacesWriteNameStr);
  
  m_interfaceNameStr = vector<string>();
  setParameter("InterfacesNames",&m_interfaceNameStr);
  
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
  
  m_fullConfigure = 1;
  
  cf_assert (m_coupledNamespacesStr.size() > 0);
  
  if (m_coupledSubSystemsStr.size() == 0) {
    const string ssname = SubSystemStatusStack::getCurrentName();
    m_coupledSubSystemsStr.resize(m_coupledNamespacesStr.size(), ssname);
  }
  
  cf_assert (m_coupledNamespacesStr.size() == m_coupledSubSystemsStr.size());
  
  // Check size of configuration parameters
  // cf_assert(m_preProcessWriteStr.size()   == m_preProcessWriteNameStr.size());
  // cf_assert(m_preProcessReadStr.size()    == m_preProcessReadNameStr.size());
  // cf_assert(m_matchMeshesWriteStr.size()  == m_matchMeshesWriteNameStr.size());
  // cf_assert(m_matchMeshesReadStr.size()   == m_matchMeshesReadNameStr.size());
  // cf_assert(m_interfacesReadStr.size()    == m_interfacesReadNameStr.size());
  cf_assert(m_interfacesWriteStr.size()   == m_interfacesWriteNameStr.size());
  // cf_assert(m_postProcessStr.size()       == m_postProcessNameStr.size());
  // cf_assert(m_coupledSubSystemsStr.size() == m_interfaceNameStr.size());

  // //Set default value for m_coupledNamespacesStr
  // if (m_coupledNamespacesStr.size() == 0) {
  //   m_coupledNamespacesStr.resize(m_coupledSubSystemsStr.size(), "Default");
  // }
  // cf_assert(m_coupledNamespacesStr.size() == m_interfaceNameStr.size());
  
  // // Set default value for m_transferRates
  // if(m_transferRates.size() == 0) {
  //   m_transferRates.resize(m_coupledSubSystemsStr.size(), 1);
  // }
  
  // cf_assert(m_transferRates.size() == m_interfaceNameStr.size());
  
  // cf_assert(m_preProcessWriteStr.size()   == m_groups.size());
  // cf_assert(m_preProcessReadStr.size()    == m_groups.size());
  // cf_assert(m_matchMeshesWriteStr.size()  == m_groups.size());
  // cf_assert(m_matchMeshesReadStr.size()   == m_groups.size());
  // cf_assert(m_interfacesReadStr.size()    == m_groups.size());
  // cf_assert(m_interfacesWriteStr.size()   == m_groups.size());
  // cf_assert(m_postProcessStr.size()       == m_groups.size());
  
  // cf_assert(m_coupledSubSystemsStr.size() == m_groups.size());
  // cf_assert(m_coupledNamespacesStr.size() == m_groups.size());
  // cf_assert(m_transferRates.size()        == m_groups.size());

  configureNested ( m_data.getPtr(), args );

  // configureInterfaces();

  // add configures to the ConcurrentCouplerCom's
  clearComds();

  // Make the name accessible to all commands by putting them in the data
  // m_data->setCoupledSubSystemsNames(m_coupledSubSystemsStr, m_coupledNamespacesStr);
  // m_data->resizeCoupledInterfaces();
  
  // // Resize the vectors of commands
  // m_preProcessRead.resize(m_preProcessReadStr.size());
  // m_preProcessWrite.resize(m_preProcessWriteStr.size());
  // m_matchMeshesRead.resize(m_matchMeshesReadStr.size());
  // m_matchMeshesWrite.resize(m_matchMeshesWriteStr.size());
  // m_interfacesRead.resize(m_interfacesReadStr.size());
  m_interfacesWrite.resize(m_interfacesWriteStr.size());
  // m_postProcess.resize(m_postProcessStr.size());
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::postConfigure( Config::ConfigArgs& args )
{
  m_fullConfigure++;
  
  // m_data->resizeCoupledInterfacesTRS();
  
  configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
    (args, m_setup,m_setupStr,m_setupStr,m_data);
  
  // Configure all commands
  for(CFuint i = 0; i < m_groups.size(); ++i) {
  //   CFLog(INFO, "ConcurrentCoupler [" << getName() << "] Namespace [" << getNamespace()
  // 	  << "] Data namespace [" << getMethodData()->getNamespace() << "]\n");
  //   CFLog(INFO, "Interface name = "    << m_interfaceNameStr[i] << "\n");
  //   CFLog(INFO, "PreProcess Reader type = "   << m_preProcessReadStr[i]    << "\n");
  //   CFLog(INFO, "PreProcess Writer type = "   << m_preProcessWriteStr[i]    << "\n");
  //   CFLog(INFO, "Mesh Matcher Reader type = " << m_matchMeshesReadStr[i]   << "\n");
  //   CFLog(INFO, "Mesh Matcher Writer type = " << m_matchMeshesWriteStr[i]   << "\n");
  //   CFLog(INFO, "Interface Reader type = " << m_interfacesReadStr[i] << "\n");
  CFLog(INFO, "Interface Writer type = " << m_interfacesWriteStr[i] << "\n");
  //   CFLog(INFO, "PostProcess type = "  << m_postProcessStr[i]   << "\n");
    
  //   CFLog(INFO, "Configuring PreProcess Reader command ["  << m_preProcessReadStr[i]   << "]\n");
  //   configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
  //     (args,m_preProcessRead[i],m_preProcessReadStr[i],m_preProcessReadNameStr[i],m_data);
    
  //   CFLog(INFO, "Configuring PreProcess Writer command ["  << m_preProcessWriteStr[i]   << "]\n");
  //   configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
  //     (args,m_preProcessWrite[i],m_preProcessWriteStr[i],m_preProcessWriteNameStr[i],m_data);
    
  //   CFLog(INFO, "Configuring Mesh Matcher Reader command ["  << m_matchMeshesReadStr[i]   << "]\n");
  //   configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
  //     (args,m_matchMeshesRead[i],m_matchMeshesReadStr[i],m_matchMeshesReadNameStr[i],m_data);
    
  //   CFLog(INFO, "Configuring Mesh Matcher Writer command ["  << m_matchMeshesWriteStr[i]   << "]\n");
  //   configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
  //     (args,m_matchMeshesWrite[i],m_matchMeshesWriteStr[i],m_matchMeshesWriteNameStr[i],m_data);
    
  //   CFLog(INFO, "Configuring Interface Reader command ["  << m_interfacesReadStr[i]   << "]\n");
  //   configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
  //     (args,m_interfacesRead[i],m_interfacesReadStr[i], m_interfacesReadNameStr[i], m_data);
    
  CFLog(INFO, "Configuring Interface Writer command ["  << m_interfacesWriteStr[i]   << "]\n");
  configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
    (args,m_interfacesWrite[i],m_interfacesWriteStr[i], m_interfacesWriteNameStr[i], m_data);
    
  //   CFLog(INFO, "Configuring PostProcess command ["  << m_postProcessNameStr[i]   << "]\n");
  //   configureCommand<ConcurrentCouplerCom,ConcurrentCouplerData,ConcurrentCouplerComProvider>
  //     (args,m_postProcess[i],m_postProcessStr[i],m_postProcessNameStr[i],m_data);
    
  //   // Check that configuration was ok
  //   cf_assert(m_preProcessRead[i].isNotNull());
  //   cf_assert(m_preProcessWrite[i].isNotNull());

  //   cf_assert(m_matchMeshesRead[i].isNotNull());
  //   cf_assert(m_matchMeshesWrite[i].isNotNull());

  //   cf_assert(m_interfacesRead[i].isNotNull());
  cf_assert(m_interfacesWrite[i].isNotNull());

  //   cf_assert(m_postProcess[i].isNotNull());
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

  m_fullConfigure++;
  
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

  m_fullConfigure++;
  
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
  //  cf_assert(m_fullConfigure == 4);
  
  cf_assert(m_setup.isNotNull());
  m_setup->execute();
  
  setupCommandsAndStrategies();
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

  cf_assert(isSetup());
  cf_assert(isConfigured());

  // const string nameSpace = getMethodData()->getNamespace();
 //  const bool isParallel = Common::PE::GetPE().IsParallel();
 //  if(isParallel)
 // {
 //   Common::PE::GetPE().setBarrier(nameSpace);
 // }

 //  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
 //  for(CFuint i = 0; i < m_interfacesRead.size(); ++i) {
 //    cf_assert(m_interfacesRead[i].isNotNull());

 //    // Execute and save file if needed...
 //    if((!(iter % m_transferRates[i])) || (iter ==0) ) {
 //      m_interfacesRead[i]->execute();
 //    }

 //  }
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::dataTransferWriteImpl()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "ConcurrentCouplerMethod::dataTransferWriteImpl() => start\n");
  
  cf_assert(isSetup());
  cf_assert(isConfigured());
  
  // if () { 
  /*const string groupName = getNamespace(); // SubSystemStatusStack::getCurrentName();
  PE::GetPE().setBarrier(groupName);
  cout << "ConcurrentCouplerMethod::dataTransferWriteImpl() => P" << PE::GetPE().GetRank(groupName) << "\n";
  
  for (;;) {}*/
  // }
  
  // const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // ///@todo this should be done only if needed
  // getMethodData()->getCollaborator<SpaceMethod>()->extrapolateStatesToNodes();
  
  for(CFuint i = 0; i < m_interfacesWrite.size(); ++i) {

    cf_assert(m_interfacesWrite[i].isNotNull());
    
    //   // Execute and save file if needed...
    //   if((!(iter % m_transferRates[i])) || (iter ==0) ) {
    m_interfacesWrite[i]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::unsetMethodImpl()
{
  // for(CFuint i = 0; i < m_postProcess.size(); ++i)
  // {
  //   cf_assert(m_postProcess[i].isNotNull());
  //   m_postProcess[i]->execute();
  // }

  unsetupCommandsAndStrategies();

  CouplerMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerMethod::clearComds()
{
  CFAUTOTRACE;

  vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_preProcessRead);
  vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_preProcessWrite);
  vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_matchMeshesRead);
  vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_matchMeshesWrite);
  vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_interfacesRead);
  vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_interfacesWrite);
  vector<SelfRegistPtr<ConcurrentCouplerCom> >().swap(m_postProcess);

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

    } // namespace ConcurrentCoupler 

  } // namespace Numerics

} // namespace COOLFluiD

