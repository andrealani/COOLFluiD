#include <QtXml>
#include <QHostInfo>
#include <mpi.h>

#include "Common/xmlParser.h"
#include "Common/PE.hh"

#include "Config/ConfigArgs.hh"
#include "Config/ConverterTools.hh"
#include "Config/ConfigFileReader.hh"

#include "ClientServer/server/ServerNetworkComm.hh"
#include "ClientServer/server/ServerSimulation.hh"
#include "ClientServer/server/SimulationManager.hh"
#include "ClientServer/server/TypesNotFoundException.hh"

#include "ClientServer/server/ServerKernel.hh"

using namespace MPI;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::network;
using namespace COOLFluiD::server;

ServerKernel::ServerKernel(const QString & hostname, quint16 portNumber, 
                           const QList<HostInfos> & hostList)
: DEFAULT_PATH(".")
{
  m_commServer = new ServerNetworkComm(hostname, portNumber);
  
  m_fileOpen = false;
  m_simRunning = false;
  m_active = false;
  
  m_hostList = hostList;
  
  // load available types from XML file
  QFile typesFile("./TypesList.xml"); // list of available types
  
  if (typesFile.open(QIODevice::ReadOnly) && m_types.setContent(&typesFile))
    typesFile.close();
  
  else 
  {
    QString message = QString("Unable to load type list from '%1': %2")
    .arg(typesFile.fileName())
    .arg(typesFile.errorString());
    
    throw TypesNotFoundException(FromHere(), message.toStdString());
  }
  
  this->createSimulator("Simulator");
  
  qRegisterMetaType<XMLNode>("XMLNode");
  
  // connect useful signals to slots
  connect(&m_simulationManager, SIGNAL(message(const QString &)), 
          this, SLOT(message(const QString &)));
  
  connect(&m_simulationManager, SIGNAL(simulationStatus(const QString &, int, 
                                                              const QString &)), this, SLOT(simulationStatus(const QString &, int, 
                                                                                                             const QString &)));
  
  connect(&m_simulationManager, SIGNAL(simulationTree(const XMLNode &)), 
          this, SLOT(simulationTree(const XMLNode &)));
  
  connect(&m_simulationManager, SIGNAL(simulationFinished()), 
          this, SLOT(simulationFinished()));
  
  connect(&m_simulationManager, SIGNAL(spawned()), this, SLOT(spawned()));
  
  connSig(SIGNAL(newClient(int)), 
          SLOT(newClient(int)));
  
  connSig(SIGNAL(getTree(int)), 
          SLOT(getTree(int)));
  
  connSig(SIGNAL(configureSimulator(int, const QDomDocument&)), 
          SLOT(configureSimulator(int, const QDomDocument&)));
  
  connSig(SIGNAL(modifyNode(int, const QDomDocument&)), 
          SLOT(modifyNode(int, const QDomDocument&)));
  
  connSig(SIGNAL(deleteNode(int, const QString&)), 
          SLOT(deleteNode(int, const QString&)));
  
  connSig(SIGNAL(addNode(int, const QString&, const QString&, const QString&, const QString&)), 
          SLOT(addNode(int, const QString&, const QString&, const QString&, const QString&)));
  
  connSig(SIGNAL(renameNode(int, const QString&, const QString&)), 
          SLOT(renameNode(int, const QString&, const QString&)));
  
  connSig(SIGNAL(getAbstractTypes(int, const QString&)), 
          SLOT(getAbstractTypes(int, const QString&)));
  
  connSig(SIGNAL(getConcreteTypes(int, const QString&)), 
          SLOT(getConcreteTypes(int, const QString&)));
  
  connSig(SIGNAL(createDirectory(int, const QString&, const QString&)), 
          SLOT(createDirectory(int, const QString&, const QString&)));
  
  connSig(SIGNAL(openDirectory(int, const QString&, const QStringList&, bool, bool)),
          SLOT(openDirectory(int, const QString&, const QStringList&, bool, bool)));
  
  connSig(SIGNAL(saveConfiguration(int, const QString&, const QDomDocument&)), 
          SLOT(saveConfiguration(int, const QString&, const QDomDocument&)));
  
  connSig(SIGNAL(openFile(int, const QString&)), 
          SLOT(openFile(int, const QString&)));
  
  connSig(SIGNAL(closeFile(int)), 
          SLOT(closeFile(int)));
  
  connSig(SIGNAL(runSimulation(int)), 
          SLOT(runSimulation(int)));
  
  connSig(SIGNAL(activateSimulation(int, unsigned int, const QString &)), 
          SLOT(activateSimulation(int, unsigned int, const QString &)));
  
  connSig(SIGNAL(deactivateSimulation(int)), 
          SLOT(deactivateSimulation(int)));
  
  connSig(SIGNAL(shutdownServer(int)), 
          SLOT(shutdownServer(int)));
  
  connSig(SIGNAL(getHostList(int)), 
          SLOT(getHostList(int))); 
  
  connSig(SIGNAL(getSubSysList(int)), 
          SLOT(getSubSysList(int)));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ServerKernel::~ServerKernel()
{
  delete m_commServer;
  delete m_srvSimulation;
}

/****************************************************************************
 
 PRIVATE METHODS
 
 ****************************************************************************/

void ServerKernel::createSimulator(const QString & name)
{
  m_srvSimulation = new ServerSimulation(name);
  
  connect(m_srvSimulation, SIGNAL(message(const QString&)), this,
          SLOT(message(const QString&)));
  
  connect(m_srvSimulation, SIGNAL(error(const QString&)), this,
          SLOT(error(const QString&)));
  
  connect(m_srvSimulation, SIGNAL(m_finished()),
          this, SLOT(simulationFinished()));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerKernel::connSig(const char * signal, const char * slot)
{
  return connect(m_commServer, signal, this, slot);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::setStatus(WorkerStatus::Type status)
{
  PE::setCurrentStatus(status);
  //  this->commServer->sendStatus(-1, PE::getStatusString(status).c_str());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerKernel::getDirContent(const QString & directory,
                                 const QStringList & extensions,
                                 bool includeFiles,
                                 bool includeNoExtension,
                                 QStringList & dirsList, 
                                 QStringList & filesList) const
{
  QStringList list;
  QDir dir(directory);
  
  dir.setFilter(QDir::Files | QDir::Dirs | QDir::Hidden);
  dir.setSorting(QDir::DirsFirst | QDir::Name);
  
  if(!dir.exists())
    return false;
  
  QFileInfoList files = dir.entryInfoList();
  QFileInfoList::iterator it = files.begin();
  
  QRegExp regex("", Qt::CaseSensitive, QRegExp::RegExp);
  
  if(!extensions.isEmpty())
  {
    /* build the regex pattern string. 
     For example, if the QStringList contains "xml" and "CFcase" extensions, 
     the resulting string will be : "^.+\\.((xml)|(CFcase))$" */
    
    QString regexPattern = extensions.join(")|(");
    regexPattern.prepend("^.+\\.((").append("))$");
    regex.setPattern(regexPattern);
  }
  else
    regex.setPattern("^.+\\..+$");
  
  while(it != files.end())
  {
    QFileInfo fileInfo = *it;
    QString filename = fileInfo.fileName();
    
    if (filename != "." && filename != "..")
    {
      if(fileInfo.isDir())
        dirsList << filename;
      
      else if(includeFiles && regex.exactMatch(filename))
        filesList << filename;
      
      else if(includeFiles && includeNoExtension && !filename.contains('.'))
        filesList << filename;
    }
    it++;
  }
  
  return true;
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void ServerKernel::newClient(int clientId)
{
  if(m_fileOpen)
    m_commServer->sendAck(clientId, true, NETWORK_OPEN_FILE);
  
  if(m_simRunning)
    m_commServer->sendAck(clientId, true, NETWORK_SIMULATION_RUNNING);
  
  if(m_active)
    m_commServer->sendAck(clientId, true, NETWORK_ACTIVATE_SIMULATION);
  
  // send a welcome message to the new client
  m_commServer->sendMessage(clientId, "Welcome to the Client-Server project!");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::getTree(int clientId)
{
  QDomDocument d; 
  d.setContent(m_srvSimulation->getTreeXML());
  m_commServer->sendTree(clientId, d);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::configureSimulator(int clientId, const QDomDocument & config)
{
  m_simTree = ConverterTools::xmlToXCFcase(config.toString().toStdString()); 
  //  ConfigArgs args = ConverterTools::xmlToConfigArgs(config.toString().toStdString()); 
  QString tree = m_simTree.createXMLString();
  QDomDocument document;
  
  document.setContent(tree);
  
  if(m_fileOpen)
    this->closeFile(clientId);
  
  m_srvSimulation->configureSimulator(document);
  
  m_fileOpen = true;
  m_commServer->sendAck(clientId, true, NETWORK_OPEN_FILE);
  
  // send the new tree to all client
  document.setContent(m_srvSimulation->getTreeXML());
  m_commServer->sendTree(-1, document);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::modifyNode(int clientId, const QDomDocument & data)
{
  if(!m_fileOpen)
    m_commServer->sendError(clientId, "No case file loaded !");
  
  else
    ; /// @todo forward to the simulator
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::deleteNode(int clientId, const QString & nodePath)
{
  if(!m_fileOpen)
    m_commServer->sendError(clientId, "No case file loaded !");
  
  else if(m_simRunning)
    m_commServer->sendError(clientId, "A simulation is running.");\
  
  else
    ; /// @todo forward to the simulator
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::addNode(int clientId, const QString & parentPath, 
                           const QString & name, const QString & type, 
                           const QString & absType)
{
  if(!m_fileOpen)
    m_commServer->sendError(clientId, "No case file loaded !");
  
  else if(m_simRunning)
    m_commServer->sendError(clientId, "A simulation is running.");
  
  else
    ; /// @todo forward to the simulator
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::renameNode(int clientId, const QString & nodePath, 
                              const QString & newName)
{
  if(!m_fileOpen)
    m_commServer->sendError(clientId, "No case file loaded !");
  
  else if(m_simRunning)
    m_commServer->sendError(clientId, "A simulation is running.");
  
  else
    ; /// @todo forward to the simulator
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::getAbstractTypes(int clientId, const QString & typeName)
{
  QDomNode node = m_types.namedItem(typeName);
  QDomDocument document;
  QDomNodeList childNodes;
  QStringList typeList;
  
  // if the node is null, typeName is not a existing type man
  if(node.isNull())
  {
    m_commServer->sendError(clientId, QString("Type '%1' does not exist.")
                                  .arg(typeName));
    return;
  }
  
  childNodes = node.childNodes();
  
  // if no child, no types to send
  if(childNodes.isEmpty())
  {
    m_commServer->sendError(clientId, QString("No abstract type for type '%1'")
                                  .arg(typeName));
    return;
  }
  
  // building the types list
  for(int i = 0 ; i < childNodes.count() ; i++)
    typeList << childNodes.at(i).nodeName();
  
  m_commServer->sendAbstractTypes(clientId, typeName, typeList);
  /*
   document = FrameBuilder::buildTypesList(TYPE_ABSTRACT_TYPES, typeName, typesList);
   
   this->commServer->send(clientId, document.toString());*/
  
  // remember that this client works on this node
  //  this->clients[client] = node;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::getConcreteTypes(int clientId, const QString & typeName)
{
  QStringList typeList = m_srvSimulation->getConcreteTypes(typeName);
  m_commServer->sendConcreteTypes(clientId, typeName, typeList);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::createDirectory(int clientId, const QString & dirPath, 
                                   const QString & name)
{
  /// @todo check if absolute
  /// @todo check if OS compatible
  QDir dir(dirPath + "/" + name);
  QString message = QString("Could not create '%1': ").arg(dir.absolutePath());
  
  // if a file or directory with that name already exists
  if(dir.exists())
  {
    message.append("a file or a directory with that name already exists");
    m_commServer->sendError(clientId, message);
  }
  
  // if directory creation failed
  else if(!dir.mkpath(dir.absolutePath()))
  {
    message.append("please check that you have the permission to create a "
                   "directory there.");
    m_commServer->sendError(clientId, message);
  }
  
  else
    m_commServer->sendAck(clientId, true, NETWORK_CREATE_DIR);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::openDirectory(int clientId, const QString & dirPath, 
                                 const QStringList & extensions, 
                                 bool includeFiles, bool includeNoExtension)
{
  QStringList directories;
  QStringList files;
  bool dotDot;
  
  QString directory;
  
  QDomDocument filesList;
  
  if(dirPath.isEmpty())
    directory = this->DEFAULT_PATH;
  else
    directory = dirPath;
  
  directory = QDir(directory).absolutePath();
  directory = QDir::cleanPath(directory);
  
  if(directory != "/") 
    directories << "..";
  
  if(!this->getDirContent(directory, extensions, includeFiles, 
                          includeNoExtension, directories, files))
  {
    m_commServer->sendError(clientId, QString("'%1' is not an existing directory")
                                  .arg(directory));
    return;
  }
  
  m_commServer->sendDirContents(clientId, directory, directories, files);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::saveConfiguration(int clientId, const QString & filename, 
                                     const QDomDocument & config)
{
  bool cfcase = false;
  QFile file(filename);
  QString configStr = config.toString();
  
  if(filename.endsWith(".CFcase"))
    cfcase = true;
  
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  {
    QString error = QString("Could open file '%1' for write access: %2")
    .arg(filename)
    .arg(file.errorString());
    
    m_commServer->sendError(clientId, error);
    m_commServer->sendAck(clientId, false, NETWORK_SAVE_CONFIG);
  }
  
  else
  {
    QTextStream out;
    
    out.setDevice(&file);
    
    if(cfcase)
    {
      std::string data = ConverterTools::xmlToCFcase(configStr.toStdString());
      out << data.c_str();
    }
    
    else
      out << configStr;
    
    file.close();
    m_commServer->sendAck(clientId, true, NETWORK_SAVE_CONFIG);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::openFile(int clientId, const QString & filename)
{
  if(m_fileOpen)
    this->closeFile(clientId);
  
  QDomDocument document;
  boost::filesystem::path file(filename.toStdString());
  
  ConfigFileReader configFile;
  ConfigArgs args;
  configFile.parse(filename.toStdString(), args);
  
  m_simTree = XMLNode::parseString(ConverterTools::configArgsToXml(args).c_str());
  
  if(m_srvSimulation->loadCaseFile(filename))
  {
    // Notify all clients that a case file has been loaded
    m_commServer->sendAck(-1, true, NETWORK_OPEN_FILE);
    m_fileOpen = true;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::closeFile(int clientId)
{
  QString error;
  
  if(m_simRunning)
    error = "A simulation is running. You cannot close the file for now";
  else if(!m_fileOpen)
    error = "No file open."; 
  if(!error.isEmpty())
    m_commServer->sendError(clientId, error);
  else
  {
    delete m_srvSimulation;
    this->createSimulator("Simulator");
    m_commServer->sendMessage(-1, "File closed");
    m_commServer->sendAck(-1, true, NETWORK_CLOSE_FILE);
    m_fileOpen = false;
  }
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::runSimulation(int clientId)
{
  QString error;
  
  if(!m_fileOpen)
    error = "Please open a case file or configure the simulator before running "
    "a simulation.";
  else if(m_simRunning)
    error = "The simulation is already running. You cannot run it twice at the "
    "same time.";
  else if(!m_active)
    error = "The simulation is not active.";
  
  if(!error.isEmpty())
    m_commServer->sendError(clientId, error);
  else
  {
    m_simulationManager.configure(m_simTree);
    
    //   this->simulationManager.start();
    
    m_simRunning = true; 
    
    // Notify all clients that the simulation has started
    m_commServer->sendAck(-1, true, NETWORK_SIMULATION_RUNNING);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::shutdownServer(int clientId)
{
  qApp->exit(0); // exit the Qt event loop
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::message(const QString & message)
{
  m_commServer->sendMessage(-1, message);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::error(const QString & message)
{
  m_commServer->sendError(-1, message);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::simulationFinished()
{
  m_simRunning = false;
  m_commServer->sendAck(-1, true, NETWORK_RUN_SIMULATION);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::getHostList(int clientId)
{
  m_commServer->sendHostList(clientId, m_hostList);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::activateSimulation(int clientId, unsigned int nbProcs, 
                                      const QString & hosts)
{
  QString error;
  
  if(!m_fileOpen)
    error = "Please open a case file or configure the simulator before activating "
    "the simulation.";
  else if(m_simRunning)
    error = "The simulation is already running.";
  else if(m_active)
    error = "The simulation is active.";
  
  if(!error.isEmpty())
  {
    m_commServer->sendError(clientId, error);
    m_commServer->sendAck(clientId, false, NETWORK_ACTIVATE_SIMULATION);
  }
  else
  {
    this->setStatus(WorkerStatus::STARTING);
    m_simulationManager.spawn("SubSystem", "SubSystem", nbProcs, hosts);
    m_active = true;
    m_simulationManager.run();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::deactivateSimulation(int clientId)
{
  QString error;
  
  if(!m_active)
    m_commServer->sendError(clientId, "The simulation is not active.");
  else
  {
    this->setStatus(WorkerStatus::EXITING);
    m_commServer->sendMessage(-1, "Exiting workers.");
    m_simulationManager.exitWorkers();
    m_commServer->sendAck(-1, true, NETWORK_DEACTIVATE_SIMULATION);
    CFout << "Simulation has been deativated.\n" << CFendl;
    m_active = false;
    this->setStatus(WorkerStatus::NOT_RUNNING);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::simulationStatus(const QString & subSysName, int rank, 
                                    const QString & status)
{
  m_commServer->sendStatus(-1, subSysName, rank, status);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::getSubSysList(int clientId)
{
  QStringList list;
  
  int size = m_srvSimulation->readSubSystems();
  
  for(int i = 0 ; i < size ; i++)
    list << m_srvSimulation->getSubSystem(i).split(" ").at(0);
  
  m_commServer->sendSubSystemList(clientId, list);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::simulationTree(const XMLNode & tree)
{
  m_commServer->sendTree(-1, tree);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerKernel::spawned()
{
  m_commServer->sendAck(-1, true, NETWORK_ACTIVATE_SIMULATION);
  m_commServer->sendMessage(-1, "Simulation has been activated.");
}