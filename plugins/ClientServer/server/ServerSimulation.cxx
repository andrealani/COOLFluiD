#include <exception>
#include <vector>
#include <string>

#include <mpi.h>

#include <QtXml>
#include <QStringList>

#include <boost/algorithm/string.hpp>

#include "logcpp/PatternLayout.hh"

#include "Common/CFLog.hh"
#include "Common/EventHandler.hh"
#include "Common/SharedPtr.hh"

#include "Config/ConverterTools.hh"
#include "Config/ConfigFileReader.hh"

#include "Environment/CFEnv.hh"
#include "Environment/FactoryRegistry.hh"
#include "Environment/FactoryBase.hh"

#include "Framework/GlobalStopCriteria.hh"
#include "Framework/Simulator.hh"
#include "Framework/SimulationStatus.hh"


#include "ClientServer/server/RemoteClientAppender.hh"
#include "ClientServer/server/ServerSimulation.hh"

using namespace MPI;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::server;
using namespace std;

ServerSimulation::ServerSimulation(const QString & simulatorName)
: Maestro(simulatorName.toStdString()), QThread()
{
  logcpp::PatternLayout * f_layout;
  logcpp::Appender * remote_appender;
  
  
  m_simulator = /*NULL;//*/new Simulator(simulatorName.toStdString());
  
  //  this->createSimulator();
  
  f_layout = new logcpp::PatternLayout();
  remote_appender = new RemoteClientAppender("RemoteClientAppender" );
  
  f_layout->setConversionPattern( "%m" );//%p 
  remote_appender->setLayout(f_layout);
  
  CFLogger::getInstance().getMainLogger().addAppender(remote_appender);
  
  m_configured = false;
  m_lastSubsystemConfigured = -1;
  
  connect((RemoteClientAppender*)remote_appender, SIGNAL(newData(const QString &)), 
          this, SLOT(newData(const QString &)));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ServerSimulation::~ServerSimulation()
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerSimulation::run()
{
  if(m_lastSubsystemConfigured == -1)
    emit error("No subsystem configured");
  
  else
  {
    try
    {   
      emit message("Starting the simulation");
      
      std::string subSysName;
      std::string subSysType;
      std::string subSys;
      std::vector<std::string> subsysnames;
      QStringList::iterator it = m_subsystemNames.begin();
      Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
      
      SimulationStatus& simStatus = SimulationStatus::getInstance();
      
      subSysName = m_subsystemNames.at(m_lastSubsystemConfigured).toStdString();
      subSysType = m_subsystemTypes.at(m_lastSubsystemConfigured).toStdString();
      subSys = subSysName + '\n' + subSysType + '\n';
      
      while(it != m_subsystemNames.end())
      {
        subsysnames.push_back(it->toStdString());
        it++;
      }
      
      simStatus.resetAll();
      simStatus.setSubSystems(subsysnames);
      
      if(m_stopcriteria.isNull())
      {
        emit error("Stop criteria is NULL...");
        return;
      }
      
      COMM_WORLD.Barrier();
      CFout << "#\n###### RUN PHASE ####################\n#\n";
      for ( ; !m_stopcriteria->isSatisfied(); )
      {
        simStatus.incrementNbIter();
        event_handler->call_signal ( "CF_ON_MAESTRO_RUN", subSys );
      }
      
      emit message("Simulation finished");
    }
    catch ( std::exception& e )
    {
      emit error(e.what());
    }
    catch (...)
    {
      emit error("Unknown exception thrown and not caught !!!\nAborting ...");
    }
    
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerSimulation::loadCaseFile(const QString & filename)
{
  try
  {
    ConfigArgs args;
    ConfigFileReader reader;
    
    reader.parse(filename.toStdString(), args);
    
    //   this->createSimulator();
    
    m_simulator->openCaseFile(filename.toStdString());
    this->configure(args);
    m_caseFile = filename;
    
    emit message("File loaded: " + filename);
    
    return true;
  }
  catch ( std::exception& e )
  {
    emit error(e.what());
  }
  catch (...)
  {
    emit error("Unknown exception thrown and not caught !!!\nAborting ...");
  }
  return false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerSimulation::newData(const QString & data)
{
  emit message(data);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QStringList ServerSimulation::getConcreteTypes(const QString & abstractType)
const
{
  QStringList typesList;
  
  /// @bug what if the abstract type does not exist ????
  
  std::vector< ProviderBase* > registered_providers =
  CFEnv::getInstance().getFactoryRegistry()->
  getFactory(abstractType.toStdString())->getAllProviders();
  
  for(size_t i = 0; i < registered_providers.size(); ++i)
    typesList << QString(registered_providers[i]->getProviderName().c_str());
  
  return typesList;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerSimulation::configureSimulator(const QDomDocument & document)
{
  QString filename;
  QTemporaryFile tempFile;
  QTextStream out;
  std::string tree;
  ConfigArgs args;
  QDomDocument doc = document;
  
  //  this->createSimulator();
  
  if(doc.firstChild().nodeName() != "XCFcase")
  {
    QDomElement elem = doc.createElement("XCFcase");
    
    elem.appendChild(doc.firstChild());
    doc.appendChild(elem);
  }
  
  tree = doc.toString().toStdString();
  args = ConverterTools::xmlToConfigArgs(tree);
  
  tempFile.open();
  filename = tempFile.fileName();
  
  out.setDevice(&tempFile);
  
  out << ConverterTools::xmlToCFcase(tree).c_str();
  
  tempFile.close();
  
  this->loadCaseFile(filename);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString ServerSimulation::getTreeXML() const
{
  return m_simulator->getTreeXML().c_str();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int ServerSimulation::getSubsystemsCount() const
{
  return m_subsystemNames.size();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerSimulation::runConfigPhase(int subsystem)
{
  bool configured = false;
  
  if(subsystem > m_subsystemNames.size())
    emit error("Unknown subsystem"); 
  else
  {
    try
    {
      std::string subSysName = m_subsystemNames.at(subsystem).toStdString();
      std::string subSysType = m_subsystemTypes.at(subsystem).toStdString();
      std::string subSys = subSysName + '\n' + subSysType + '\n';
      Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
      
      // we manage the simulator
      this->manage(m_simulator);
      
      CFout << "+++++++++++++++++++++++++++++++ CONFIGURING with " << __FUNCTION__
      << " ++++++++++++++++++++++++++++++++\n" << CFendl;
      
      CFout << "#\n###### STARTING SUBSYSTEM : " << subSysName << " ######\n#\n";
      CFout << "in " << __FUNCTION__ << " at line " << __LINE__ << "\n" << CFendl;
      event_handler->call_signal("CF_ON_MAESTRO_BUILDSUBSYSTEM", subSys);
      
      //    COMM_WORLD.Barrier();
      
      CFout << "#\n###### CONFIG PHASE #################\n#\n";
      event_handler->call_signal("CF_ON_MAESTRO_CONFIGSUBSYSTEM", subSys);
      
      //    COMM_WORLD.Barrier();
      
      CFout << "#\n###### SOCKETS PLUG PHASE ###########\n#\n";
      event_handler->call_signal("CF_ON_MAESTRO_PLUGSOCKETS", subSys);
      
      //    COMM_WORLD.Barrier();
      
      CFout << "#\n###### BUILD PHASE ##################\n#\n";
      event_handler->call_signal("CF_ON_MAESTRO_BUILDMESHDATA", subSys);
      
      //    COMM_WORLD.Barrier();
      
      CFout << "#\n###### SETUP PHASE ##################\n#\n";
      event_handler->call_signal("CF_ON_MAESTRO_SETUP", subSys);
      
      COMM_WORLD.Barrier();
      
      m_lastSubsystemConfigured = subsystem;
      
      configured = true;
    }
    catch ( std::exception& e )
    {
      emit error(e.what());
    }
    catch (...)
    {
      emit error("Unknown exception thrown and not caught !!!\nAborting ...");
    }
  }
  return configured;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerSimulation::runUnconfigPhase()
{
  bool unconfigured = false;
  
  if(m_lastSubsystemConfigured == -1)
    emit error("No subsystem configured"); 
  else
  {
    try
    {
      std::string subSysName = m_subsystemNames.at(m_lastSubsystemConfigured).toStdString();
      std::string subSysType = m_subsystemTypes.at(m_lastSubsystemConfigured).toStdString();
      std::string subSys = subSysName + '\n' + subSysType + '\n';
      Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();
      
      CFout << "#\n###### UNSETUP PHASE ################\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_UNSETUP", subSys );
      
      //    COMM_WORLD.Barrier();
      
      CFout << "#\n###### SOCKETS UNPLUG PHASE #########\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_UNPLUGSOCKETS", subSys );
      
      //    COMM_WORLD.Barrier();
      
      CFout << "#\n###### DESTRUCTION SUBSYSTEM PHASE #########\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_DESTROYSUBSYSTEM", subSys );
      
      COMM_WORLD.Barrier();
      
      m_lastSubsystemConfigured = -1;
      
      unconfigured = true;
    }
    catch ( std::exception& e )
    {
      emit error(e.what());
    }
    catch (...)
    {
      emit error("Unknown exception thrown and not caught !!!\nAborting ...");
    }
  }
  
  return unconfigured;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int ServerSimulation::readSubSystems()
{
  // which subsystems will be controlled by meastro
  std::vector<std::string> subsysnames = m_simulator->getSubSystemNames();
  std::vector<std::string> subsystypes = m_simulator->getSubSystemTypes();
  
  cf_assert (subsysnames.size() == subsystypes.size());
  
  m_subsystemNames.clear();
  m_subsystemTypes.clear();
  
  for ( CFuint i = 0; i < subsysnames.size(); ++i )
  {
    m_subsystemNames << subsysnames[i].c_str();
    m_subsystemTypes << subsystypes[i].c_str();
    
  }
  
  //  this->subsystemNames << "Fluid";
  //  this->subsystemTypes << "Fluid";
  //  this->subsystemNames << "Mesh";
  //  this->subsystemTypes << "Mesh";
  //  this->subsystemNames << "Solid";
  //  this->subsystemTypes << "Solid";
  
  return m_subsystemNames.size();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString ServerSimulation::getSubSystem(int subSystem) const
{
  return m_subsystemNames.at(subSystem) + " " + m_subsystemTypes.at(subSystem);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerSimulation::createSimulator()
{
  //  if(this->simulator.isNotNull())
  //   this->simulator.release();
  //  
  m_simulator = new Simulator("Simulator");
}
