#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/CouplerMethod.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/DirPaths.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/Coupler.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Coupler,
               CouplerMethod,
               SubSystemCouplerModule,
               1>
couplerProvider("SubSystemCoupler");

//////////////////////////////////////////////////////////////////////////////

void Coupler::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("SetupCom","SetupCommand to run. This command seldomly needs overriding.");

   options.addConfigOption< std::vector<std::string> >("PreProcessReadComs","Preprocessing");
   options.addConfigOption< std::vector<std::string> >("PreProcessReadNames","Names for the configuration of the commands.");
   options.addConfigOption< std::vector<std::string> >("PreProcessWriteComs","Preprocessing");
   options.addConfigOption< std::vector<std::string> >("PreProcessWriteNames","Names for the configuration of the commands.");

   options.addConfigOption< std::vector<std::string> >("MeshMatchingReadComs","MeshMatching command.");
   options.addConfigOption< std::vector<std::string> >("MeshMatchingReadNames","Names for the configuration of the commands.");
   options.addConfigOption< std::vector<std::string> >("MeshMatchingWriteComs","MeshMatching command.");
   options.addConfigOption< std::vector<std::string> >("MeshMatchingWriteNames","Names for the configuration of the commands.");

   options.addConfigOption< std::vector<std::string> >("InterfacesReadNames","Names for the configuration of the commands.");
   options.addConfigOption< std::vector<std::string> >("InterfacesReadComs","dataTransreRead commands.");
   options.addConfigOption< std::vector<std::string> >("InterfacesWriteNames","Names for the configuration of the commands.");
   options.addConfigOption< std::vector<std::string> >("InterfacesWriteComs","dataTransreWrite commands.");

   options.addConfigOption< std::vector<std::string> >("PostProcessNames","Names for the configuration of the commands.");
   options.addConfigOption< std::vector<std::string> >("PostProcessComs","Unsetup commands.");

   options.addConfigOption< std::vector<std::string> >("InterfacesNames","Names of the interfaces.");
   options.addConfigOption< std::vector<std::string> >("CoupledSubSystems","Names of the subsystems to be coupled to.");
   options.addConfigOption< std::vector<std::string> >("CoupledNameSpaces","Names of the namespaces of the other subsystems.");

   options.addConfigOption< std::vector<CFuint> >("TransferRates","Transfer data every X iterations");
}

//////////////////////////////////////////////////////////////////////////////

Coupler::Coupler(const std::string& name) :
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
  m_data.reset(new SubSysCouplerData(this));

  m_setupStr = "StdSetup";
  setParameter("SetupCom",&m_setupStr);

  m_preProcessReadStr = vector<std::string>();
   setParameter("PreProcessReadComs",&m_preProcessReadStr);
  m_preProcessReadNameStr = vector<std::string>();
   setParameter("PreProcessReadNames",&m_preProcessReadNameStr);

  m_preProcessWriteStr = vector<std::string>();
   setParameter("PreProcessWriteComs",&m_preProcessWriteStr);
  m_preProcessWriteNameStr = vector<std::string>();
   setParameter("PreProcessWriteNames",&m_preProcessWriteNameStr);

  m_matchMeshesReadStr = vector<std::string>();
   setParameter("MeshMatchingReadComs",&m_matchMeshesReadStr);
  m_matchMeshesReadNameStr = vector<std::string>();
   setParameter("MeshMatchingReadNames",&m_matchMeshesReadNameStr);

  m_matchMeshesWriteStr = vector<std::string>();
   setParameter("MeshMatchingWriteComs",&m_matchMeshesWriteStr);
  m_matchMeshesWriteNameStr = vector<std::string>();
   setParameter("MeshMatchingWriteNames",&m_matchMeshesWriteNameStr);

  m_postProcessStr = vector<std::string>();
   setParameter("PostProcessComs",&m_postProcessStr);
  m_postProcessNameStr = vector<std::string>();
   setParameter("PostProcessNames",&m_postProcessNameStr);

  m_interfacesReadStr = vector<std::string>();
   setParameter("InterfacesReadComs",&m_interfacesReadStr);
  m_interfacesReadNameStr = vector<std::string>();
   setParameter("InterfacesReadNames",&m_interfacesReadNameStr);

  m_interfacesWriteStr = vector<std::string>();
   setParameter("InterfacesWriteComs",&m_interfacesWriteStr);
  m_interfacesWriteNameStr = vector<std::string>();
   setParameter("InterfacesWriteNames",&m_interfacesWriteNameStr);

  m_interfaceNameStr = vector<std::string>();
   setParameter("InterfacesNames",&m_interfaceNameStr);

  m_coupledSubSystemsStr = vector<std::string>();
   setParameter("CoupledSubSystems",&m_coupledSubSystemsStr);

  m_coupledNamespacesStr = vector<std::string>();
   setParameter("CoupledNameSpaces",&m_coupledNamespacesStr);

  m_transferRates = std::vector<CFuint>();
  setParameter("TransferRates",&m_transferRates);

}

//////////////////////////////////////////////////////////////////////////////

Coupler::~Coupler()
{
  clearComds();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> Coupler::getMethodData() const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::configure ( Config::ConfigArgs& args )
{
  Framework::CouplerMethod::configure(args);

  m_fullConfigure = 1;

  // Check size of configuration parameters
  cf_assert(m_preProcessWriteStr.size()   == m_preProcessWriteNameStr.size());
  cf_assert(m_preProcessReadStr.size()    == m_preProcessReadNameStr.size());
  cf_assert(m_matchMeshesWriteStr.size()  == m_matchMeshesWriteNameStr.size());
  cf_assert(m_matchMeshesReadStr.size()   == m_matchMeshesReadNameStr.size());
  cf_assert(m_interfacesReadStr.size()    == m_interfacesReadNameStr.size());
  cf_assert(m_interfacesWriteStr.size()   == m_interfacesWriteNameStr.size());
  cf_assert(m_postProcessStr.size()       == m_postProcessNameStr.size());
  cf_assert(m_coupledSubSystemsStr.size() == m_interfaceNameStr.size());

  //Set default value for m_coupledNamespacesStr
  if(m_coupledNamespacesStr.size() == 0)
  {
    m_coupledNamespacesStr.resize(m_coupledSubSystemsStr.size());
    for(CFuint iNsp = 0; iNsp < m_coupledNamespacesStr.size(); iNsp++)
    {
      m_coupledNamespacesStr[iNsp] = "Default";
    }
  }
  cf_assert(m_coupledNamespacesStr.size() == m_interfaceNameStr.size());

  // Set default value for m_transferRates
  if(m_transferRates.size() == 0)
  {
    m_transferRates.resize(m_coupledSubSystemsStr.size());
    for(CFuint i = 0; i < m_transferRates.size(); ++i)
    {
      m_transferRates[i] = 1;
    }
  }
  cf_assert(m_transferRates.size() == m_interfaceNameStr.size());

  cf_assert(m_preProcessWriteStr.size()   == m_groups.size());
  cf_assert(m_preProcessReadStr.size()    == m_groups.size());
  cf_assert(m_matchMeshesWriteStr.size()  == m_groups.size());
  cf_assert(m_matchMeshesReadStr.size()   == m_groups.size());
  cf_assert(m_interfacesReadStr.size()    == m_groups.size());
  cf_assert(m_interfacesWriteStr.size()   == m_groups.size());
  cf_assert(m_postProcessStr.size()       == m_groups.size());

  cf_assert(m_coupledSubSystemsStr.size() == m_groups.size());
  cf_assert(m_coupledNamespacesStr.size() == m_groups.size());
  cf_assert(m_transferRates.size()        == m_groups.size());

  configureNested ( m_data.getPtr(), args );

  configureInterfaces();

  // add configures to the CouplerCom's
  clearComds();

  // Make the name accessible to all commands by putting them in the data
  m_data->setCoupledSubSystemsNames(m_coupledSubSystemsStr, m_coupledNamespacesStr);
  m_data->resizeCoupledInterfaces();

  // Resize the vectors of commands
  m_preProcessRead.resize(m_preProcessReadStr.size());
  m_preProcessWrite.resize(m_preProcessWriteStr.size());
  m_matchMeshesRead.resize(m_matchMeshesReadStr.size());
  m_matchMeshesWrite.resize(m_matchMeshesWriteStr.size());
  m_interfacesRead.resize(m_interfacesReadStr.size());
  m_interfacesWrite.resize(m_interfacesWriteStr.size());
  m_postProcess.resize(m_postProcessStr.size());

}

//////////////////////////////////////////////////////////////////////////////

void Coupler::postConfigure( Config::ConfigArgs& args )
{

  m_fullConfigure++;
  m_data->resizeCoupledInterfacesTRS();

  configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args, m_setup,m_setupStr,m_setupStr,m_data);

  // Configure all commands
  for(CFuint i = 0; i < m_groups.size(); ++i)
  {

    CFLog(INFO, "Coupler [" << getName()
             << "] Namespace [" << getNamespace()
             << "] Data namespace [" << getMethodData()->getNamespace()
             << "]\n");

    CFLog(INFO, "Interface name = "    << m_interfaceNameStr[i] << "\n");
    CFLog(INFO, "PreProcess Reader type = "   << m_preProcessReadStr[i]    << "\n");
    CFLog(INFO, "PreProcess Writer type = "   << m_preProcessWriteStr[i]    << "\n");
    CFLog(INFO, "Mesh Matcher Reader type = " << m_matchMeshesReadStr[i]   << "\n");
    CFLog(INFO, "Mesh Matcher Writer type = " << m_matchMeshesWriteStr[i]   << "\n");
    CFLog(INFO, "Interface Reader type = " << m_interfacesReadStr[i] << "\n");
    CFLog(INFO, "Interface Writer type = " << m_interfacesWriteStr[i] << "\n");
    CFLog(INFO, "PostProcess type = "  << m_postProcessStr[i]   << "\n");


    CFLog(INFO, "Configuring PreProcess Reader command ["  << m_preProcessReadStr[i]   << "]\n");
    configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args,m_preProcessRead[i],m_preProcessReadStr[i],m_preProcessReadNameStr[i],m_data);
    CFLog(INFO, "Configuring PreProcess Writer command ["  << m_preProcessWriteStr[i]   << "]\n");
    configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args,m_preProcessWrite[i],m_preProcessWriteStr[i],m_preProcessWriteNameStr[i],m_data);
    CFLog(INFO, "Configuring Mesh Matcher Reader command ["  << m_matchMeshesReadStr[i]   << "]\n");
    configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args,m_matchMeshesRead[i],m_matchMeshesReadStr[i],m_matchMeshesReadNameStr[i],m_data);
    CFLog(INFO, "Configuring Mesh Matcher Writer command ["  << m_matchMeshesWriteStr[i]   << "]\n");
    configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args,m_matchMeshesWrite[i],m_matchMeshesWriteStr[i],m_matchMeshesWriteNameStr[i],m_data);
    CFLog(INFO, "Configuring Interface Reader command ["  << m_interfacesReadStr[i]   << "]\n");
    configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args,m_interfacesRead[i],m_interfacesReadStr[i], m_interfacesReadNameStr[i], m_data);
    CFLog(INFO, "Configuring Interface Writer command ["  << m_interfacesWriteStr[i]   << "]\n");
    configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args,m_interfacesWrite[i],m_interfacesWriteStr[i], m_interfacesWriteNameStr[i], m_data);
    CFLog(INFO, "Configuring PostProcess command ["  << m_postProcessNameStr[i]   << "]\n");
    configureCommand<CouplerCom,SubSysCouplerData,CouplerComProvider>(args,m_postProcess[i],m_postProcessStr[i],m_postProcessNameStr[i],m_data);

    // Check that configuration was ok
    cf_assert(m_preProcessRead[i].isNotNull());
    cf_assert(m_preProcessWrite[i].isNotNull());

    cf_assert(m_matchMeshesRead[i].isNotNull());
    cf_assert(m_matchMeshesWrite[i].isNotNull());

    cf_assert(m_interfacesRead[i].isNotNull());
    cf_assert(m_interfacesWrite[i].isNotNull());

    cf_assert(m_postProcess[i].isNotNull());
  }


  // Check that all the commands apply to the same group of TRSs
  for(CFuint i = 0; i < m_interfacesRead.size(); ++i)
  {
     cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_matchMeshesRead[i]->getCommandGroupName());
     cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_matchMeshesWrite[i]->getCommandGroupName());
     cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_preProcessRead[i]->getCommandGroupName());
     cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_preProcessWrite[i]->getCommandGroupName());
     cf_assert(m_interfacesRead[i]->getCommandGroupName() == m_postProcess[i]->getCommandGroupName());
     cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_matchMeshesRead[i]->getCommandGroupName());
     cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_matchMeshesWrite[i]->getCommandGroupName());
     cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_preProcessRead[i]->getCommandGroupName());
     cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_preProcessWrite[i]->getCommandGroupName());
     cf_assert(m_interfacesWrite[i]->getCommandGroupName() == m_postProcess[i]->getCommandGroupName());
  }

}

//////////////////////////////////////////////////////////////////////////////

void Coupler::configureInterfaces()
{
  // All Command Groups are assumed to be interfaces
  m_data->setInterfaces(getCommandGroups());

  // prepare the InterfacesMap vector
  SafePtr<CFMap<std::string, CFuint> > coupledInterfacesMap = m_data->getCoupledInterfacesMapPtr();
  coupledInterfacesMap->reserve(m_groups.size());

  // we need to have at least one group
  cf_assert(getNbGroups() > 0);

  // fill the map
  for(CFuint i = 0; i < m_groups.size(); ++i) {
    coupledInterfacesMap->insert(m_groups[i]->getName(),i);
  }

  coupledInterfacesMap->sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::setInfoToOtherSubSystem()
{
  CFAUTOTRACE;

  m_fullConfigure++;

  /// @todo this will be substituted by a communication stream
  const bool isParallel = Common::PE::GetPE().IsParallel();
  if((!isParallel) || ((isParallel) && (Common::PE::GetPE().GetRank() == 0)))
  {
    const std::string nameSpace = getNamespace();
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(nameSpace);
    Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
    Common::SafePtr<PhysicalModel> physicalModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
    const std::string currentSubSystem = subsystemStatus->getSubSystemName();

    // For each interface, write the number and names of TRS
    for(CFuint i = 0; i < m_groups.size(); ++i) {

      const std::string interfaceName = m_groups[i]->getName();

      const vector<std::string>& trsNames = m_groups[i]->getTrsNames();
      const CFuint nbTRS = trsNames.size();

      boost::filesystem::path dataHandleName ("COUPLING_" + interfaceName + "_" + nameSpace + "_" + currentSubSystem);
      boost::filesystem::path nameOutputFile =
        Environment::DirPaths::getInstance().getResultsDir() / dataHandleName;

      SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
      ofstream& fout = fhandle->open(nameOutputFile);

      const vector<std::string> coordType = m_data->getThisCoordType(interfaceName);
      //Get variable transformer
      Common::SafePtr<PreVariableTransformer> varTransfo = m_data->getPreVariableTransformer(interfaceName);

      fout << "!COORDTYPE " << coordType.size() << "\n";
      for (CFuint i = 0; i < coordType.size(); ++i) {
        fout << coordType[i] << "\n";
      }
      fout << "!DATASIZE " << varTransfo->getTransformedSize(physicalModel->getNbEq()) << "\n";

      fout << "!NBTRS " << nbTRS << "\n";
      for (CFuint i = 0; i < nbTRS; ++i) {
        fout << trsNames[i] << "\n";
//         CF_DEBUG_STR(trsNames[i]);
      }

      fhandle->close();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::getInfoFromOtherSubSystem()
{
  CFAUTOTRACE;
  m_fullConfigure++;

  const bool isParallel = Common::PE::GetPE().IsParallel();
  if(isParallel)
  {
    Common::PE::GetPE().setBarrier();

    for (CFuint i = 0; i < Common::PE::GetPE().GetProcessorCount(); ++i) {
      if (i == Common::PE::GetPE().GetRank()) {
        readInfoFromOtherSubSystem();
      }
      Common::PE::GetPE().setBarrier();
    }
  }
  else
  {
    readInfoFromOtherSubSystem();
  }
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::readInfoFromOtherSubSystem()
{

  // Get the information about the other subsystem TRS's coupled
  for(CFuint i = 0; i < m_groups.size(); ++i) {

    const std::string interfaceName = m_groups[i]->getName();

    const std::string otherNamespace = m_data->getCoupledNameSpaceName(interfaceName);
    const std::string otherSubSystem = m_data->getCoupledSubSystemName(interfaceName);
    boost::filesystem::path fileName ("COUPLING_" + interfaceName + "_" + otherNamespace + "_" + otherSubSystem);

    //Reading the file
    boost::filesystem::path fname =
      Environment::DirPaths::getInstance().getResultsDir() / fileName;

    Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
    ifstream& fin = fhandle->open(fname);

    CFuint lineNb = 0;
    std::string line;
    vector<std::string> words;

    // read type of Data
    getWordsFromLine(fin,line,lineNb,words);
    cf_assert(words.size() == 2);
    cf_assert(words[0] == "!COORDTYPE");
    const CFuint nbCoordType = Common::StringOps::from_str<CFint>(words[1]);
    vector<std::string> coordType(nbCoordType);
    for (CFuint i=0; i< nbCoordType; ++i)
    {
      getWordsFromLine(fin,line,lineNb,words);
      cf_assert(words.size() == 1);
      coordType[i] = words[0];
      cf_assert((coordType[i] == "Nodes")||(coordType[i] == "States") ||(coordType[i] == "Ghost") ||(coordType[i] == "Gauss"));
    }

    // read size of Data
    getWordsFromLine(fin,line,lineNb,words);
    cf_assert(words.size() == 2);
    cf_assert(words[0] == "!DATASIZE");
    const CFuint dataSize = Common::StringOps::from_str<CFint>(words[1]);

    // read nb trs's
    getWordsFromLine(fin,line,lineNb,words);
    cf_assert(words.size() == 2);
    cf_assert(words[0] == "!NBTRS");
    const CFuint nbTRS = Common::StringOps::from_str<CFint>(words[1]);
    vector<std::string> trsNames(nbTRS);

    for (CFuint i=0; i< nbTRS; ++i)
    {
      getWordsFromLine(fin,line,lineNb,words);
      cf_assert(words.size() == 1);
      trsNames[i] = words[0];
    }

    fhandle->close();

    m_data->setCoupledSubSystemTRSNames(interfaceName, trsNames);
    m_data->setCoupledSubSystemCoordTypes(interfaceName, coordType);
    m_data->setTransferedSize(interfaceName, dataSize);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::setMethodImpl()
{
  CFAUTOTRACE;

  CouplerMethod::setMethodImpl();

  //first check that the method and its commands have been correctly configured
  cf_assert(isConfigured());
  cf_assert(m_fullConfigure == 4);

  cf_assert(m_setup.isNotNull());
  m_setup->execute();

  setupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::preProcessReadImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());

  const bool isParallel = Common::PE::GetPE().IsParallel();
  if(isParallel)
  {
    Common::PE::GetPE().setBarrier();
  }

  // preprocess interfaces
  for(CFuint i = 0; i < m_preProcessRead.size(); ++i) {
    cf_assert(m_preProcessRead[i].isNotNull());
    m_preProcessRead[i]->execute();
  }

}

//////////////////////////////////////////////////////////////////////////////

void Coupler::preProcessWriteImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());

  // preprocess interfaces
  for(CFuint i = 0; i < m_preProcessWrite.size(); ++i) {
    cf_assert(m_preProcessWrite[i].isNotNull());
    m_preProcessWrite[i]->execute();
  }

}

//////////////////////////////////////////////////////////////////////////////

void Coupler::meshMatchingReadImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());

  const bool isParallel = Common::PE::GetPE().IsParallel();
  if(isParallel)
  {
    Common::PE::GetPE().setBarrier();
  }

  // match interface meshes
  for(CFuint i = 0; i < m_matchMeshesRead.size(); ++i) {
    cf_assert(m_matchMeshesRead[i].isNotNull());
    m_matchMeshesRead[i]->execute();
  }

}

//////////////////////////////////////////////////////////////////////////////

void Coupler::meshMatchingWriteImpl()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());

  // match interface meshes
  for(CFuint i = 0; i < m_matchMeshesWrite.size(); ++i) {
    cf_assert(m_matchMeshesWrite[i].isNotNull());
    m_matchMeshesWrite[i]->execute();
  }
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::dataTransferReadImpl()
{
  CFAUTOTRACE;

  cf_assert(isSetup());
  cf_assert(isConfigured());

  const bool isParallel = Common::PE::GetPE().IsParallel();
  if(isParallel)
  {
    Common::PE::GetPE().setBarrier();
  }

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  for(CFuint i = 0; i < m_interfacesRead.size(); ++i) {
    cf_assert(m_interfacesRead[i].isNotNull());

    // Execute and save file if needed...
    if((!(iter % m_transferRates[i])) || (iter ==0) ) {
      m_interfacesRead[i]->execute();
    }

  }
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::dataTransferWriteImpl()
{
  CFAUTOTRACE;

  cf_assert(isSetup());
  cf_assert(isConfigured());

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  ///@todo this should be done only if needed
  getMethodData()->getCollaborator<SpaceMethod>()->extrapolateStatesToNodes();

  for(CFuint i = 0; i < m_interfacesWrite.size(); ++i) {
    cf_assert(m_interfacesWrite[i].isNotNull());

    // Execute and save file if needed...
    if((!(iter % m_transferRates[i])) || (iter ==0) ) {
      m_interfacesWrite[i]->execute();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::unsetMethodImpl()
{
  for(CFuint i = 0; i < m_postProcess.size(); ++i)
  {
    cf_assert(m_postProcess[i].isNotNull());
    m_postProcess[i]->execute();
  }

  unsetupCommandsAndStrategies();

  CouplerMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::clearComds()
{
  CFAUTOTRACE;

  vector<SelfRegistPtr<CouplerCom> >().swap(m_preProcessRead);
  vector<SelfRegistPtr<CouplerCom> >().swap(m_preProcessWrite);
  vector<SelfRegistPtr<CouplerCom> >().swap(m_matchMeshesRead);
  vector<SelfRegistPtr<CouplerCom> >().swap(m_matchMeshesWrite);
  vector<SelfRegistPtr<CouplerCom> >().swap(m_interfacesRead);
  vector<SelfRegistPtr<CouplerCom> >().swap(m_interfacesWrite);
  vector<SelfRegistPtr<CouplerCom> >().swap(m_postProcess);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::NumericalStrategy> >
Coupler::getStrategyList() const
{
  vector<Common::SafePtr<Framework::NumericalStrategy> > result;

  std::vector < Common::SelfRegistPtr<PreVariableTransformer > > preTrans = m_data->getPreVariableTransformers();
  for(CFuint i = 0; i < preTrans.size(); ++i) {
    result.push_back(preTrans[i].getPtr());
  }

  std::vector < Common::SelfRegistPtr<PostVariableTransformer > > postTrans = m_data->getPostVariableTransformers();
  for(CFuint i = 0; i < postTrans.size(); ++i) {
    result.push_back(postTrans[i].getPtr());
  }

 return result;
}

//////////////////////////////////////////////////////////////////////////////

void Coupler::getWordsFromLine(ifstream& fin,
                               std::string& line,
                               CFuint&  lineNb,
                               vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

