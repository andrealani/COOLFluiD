
#include <QtCore>

#include "ClientServer/client/ClientNetworkComm.hh"
#include "ClientServer/client/GlobalLog.hh"
#include "ClientServer/network/NetworkFrameType.hh"
#include "ClientServer/client/StatusModel.hh"
#include "ClientServer/treeview/TSshInformation.hh"

#include "ClientServer/client/ClientKernel.hh"

#define connectSig(comm,slotSig) connect(comm, SIGNAL(slotSig), this, SLOT(slotSig));

using namespace COOLFluiD::client;
using namespace COOLFluiD::network;
using namespace COOLFluiD::treeview;

ClientKernel * ClientKernel::m_instance = NULL;

ClientKernel::ClientKernel()
{
  //  this->timer = new QTimer();
  
  m_treeModel = NULL;
  m_statusModel = NULL;
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClientKernel::~ClientKernel()
{
  //  delete this->timer;
}

/****************************************************************************
 
 PUBLIC METHODS/FUNCTIONS
 
 ****************************************************************************/

ClientKernel * ClientKernel::getInstance()
{
  if(m_instance == NULL)
    m_instance = new ClientKernel();
  
  return m_instance;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::setTreeModel(TreeModel * treeModel)
{
  if(m_treeModel != NULL)
    disconnect(m_treeModel);
  
  m_treeModel = treeModel;
  
  if(m_treeModel != NULL)
  {
    /// @todo connect signals (if any)
  }
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::setStatusModel(StatusModel * statusModel)
{
  if(m_statusModel != NULL)
    disconnect(m_statusModel);
  
  m_statusModel = statusModel;
  
  if(m_statusModel != NULL)
  {
    /// @todo connect signals (if any)
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeModel * ClientKernel::getTreeModel() const
{
  return m_treeModel;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

StatusModel * ClientKernel::getStatusModel() const
{
  return m_statusModel;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::readDir(const QModelIndex & index, const QString & dirname,
                           bool includeFiles, const QStringList & extensions,
                           bool includeNoExtension)
{
  ClientNetworkComm * comm = m_networkComms[index];
  
  if(comm == NULL)
    GlobalLog::error("Unknown simulation index");
  else
    comm->sendOpenDir(dirname, includeFiles, extensions, includeNoExtension);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::createDir(const QModelIndex & index, 
                             const QString & parentPath, const QString & dirName)
{
  ClientNetworkComm * comm = m_networkComms[index];
  
  if(comm == NULL)
    GlobalLog::error("Unknown simulation index");
  else
    comm->sendCreateDir(parentPath, dirName);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::openFile(const QModelIndex & index, const QString & filename)
{
  ClientNetworkComm * comm = m_networkComms[index];
  
  if(comm == NULL)
    GlobalLog::error("Unknown simulation index");
  else
    comm->sendOpenFile(filename);
}

/****************************************************************************
 
 PRIVATE METHODS
 
 ****************************************************************************/

QModelIndex ClientKernel::getSimIndex(QObject * sender) const
{
  ClientNetworkComm * comm = static_cast<ClientNetworkComm*>(sender);
  return m_networkComms.key(comm);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::connectToServer(const QModelIndex & simIndex)
{
  ClientNetworkComm * comm = m_networkComms[simIndex];
  TSshInformation & sshInfo = m_commsSshInfo[simIndex];
  QTimer * timer = m_timers[simIndex];
  
  if(!timer->isActive() && sshInfo.launchServer)
    this->launchServer(simIndex);
  else if(comm != NULL)
    comm->connectToServer(sshInfo.m_hostname, sshInfo.port, timer->isActive());
  else
    GlobalLog::error("Invalid index");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::launchServer(const QModelIndex & simIndex)
{
  QString cmd;
  QString msg = "Checking if no other server instance is running on "
  "\"%1\" on port \"%2\"...";
  QProcess checkIfRunning;
  TSshInformation & sshInfo = m_commsSshInfo[simIndex];
  
  //
  // checking if the remote port is available
  //
  
  cmd = QString("ssh %1@%2 check_coolfluid_server.sh %3")
  .arg(sshInfo.username)
  .arg(sshInfo.m_hostname)
  .arg(sshInfo.port);
  
  GlobalLog::message(msg.arg(sshInfo.m_hostname).arg(sshInfo.port));
  checkIfRunning.start(cmd);
  
  // wait that the process is done 
  checkIfRunning.waitForFinished(-1);
  
  QString output = checkIfRunning.readAllStandardOutput();
  QString errorMsg = checkIfRunning.readAllStandardError();
  
  // if there is something on the stderr, the check failed
  if(!errorMsg.isEmpty())
  {
    GlobalLog::error(errorMsg);
    return ;
  }
  
  // if output is different from "0", a server is already running on this port
  if(output != "0")
  {
    QString msg = "A server is already running on port %1 on %2. Please change "
    "the port or the hostname.";
    
    GlobalLog::error(msg.arg(sshInfo.m_hostname).arg(sshInfo.port));
    return ;
  }
  
  //
  // launching the server application
  //
  
  cmd = QString("ssh -n %1@%2 start_coolfluid_server.sh %3")
  .arg(sshInfo.username)
  .arg(sshInfo.m_hostname)
  .arg(sshInfo.port);
  
  GlobalLog::message("Starting the server...");
  
  // client will try every 100 ms to connect to the server
  m_timers[simIndex]->start(100);
  
  m_launchServerProcs[simIndex]->start(cmd);
}

/****************************************************************************
 
 PRIVATE SLOTS
 
 ****************************************************************************/

void ClientKernel::newTree(const QDomDocument & domDocument)
{
  QModelIndex simulation = this->getSimIndex(sender());
  
  if(m_treeModel != NULL && simulation.isValid())
  {
    m_treeModel->setSimulationTree(domDocument, simulation);  
    GlobalLog::message("Treeview updated.");
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::connected()
{
  ClientNetworkComm * comm = static_cast<ClientNetworkComm*>(sender());
  QModelIndex simulation = this->getSimIndex(this->sender());
  TSshInformation sshInfo = m_commsSshInfo[simulation];
  QString msg = "Now connected to server '%1' on port %2.";
  
  if(m_timers[simulation]->isActive())
  {
    // stop the process (send SIGKILL signal)
    m_launchServerProcs[simulation]->kill();
    m_timers[simulation]->stop();
    GlobalLog::message("Server started!");
  }
  
  m_treeModel->setSimConnectedStatus(simulation, true);
  GlobalLog::message(msg.arg(sshInfo.m_hostname).arg(sshInfo.port));
  comm->sendGetHostList();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::abstractTypes(const QStringList & types)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::concreteTypes(const QStringList & types)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::tryToConnect()
{
  QModelIndex index = m_timers.key(static_cast<QTimer*>(sender()));
  
  if(index.isValid())
  {
    TSshInformation & sshInfo = m_commsSshInfo[index];
    m_networkComms[index]->connectToServer(sshInfo.m_hostname, sshInfo.port, true);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::hostList(const QList<HostInfos> & infos)
{
  QDomDocument doc;
  QDomElement root = doc.createElement("hosts");
  QList<HostInfos>::const_iterator it = infos.begin();
  QModelIndex index = this->getSimIndex(sender());
  
  while(it != infos.end())
  {
    QDomElement elem = doc.createElement("item");
    QString text = "%1 (%2 slot(s))";
    
    elem.setAttribute("name", it->m_hostname);
    elem.setAttribute("selected", "false");
    
    elem.appendChild(doc.createTextNode(text.arg(it->m_hostname).arg(it->m_nbSlots)));
    
    root.appendChild(elem);
    
    it++;
  }
  
  doc.appendChild(root);
  
  m_treeModel->setSimulationHostList(index, doc);
  GlobalLog::message(QString("Host list updated for simulation '%1'")
                     .arg(m_treeModel->getNodePath(index)));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::subsystemList(const QStringList & subSystems)
{
  QStringList::const_iterator it = subSystems.begin();
  QModelIndex index = this->getSimIndex(sender());
  int nbProcs;
  QStringList hosts;
  
  m_treeModel->getWorkersInfo(index, nbProcs, hosts);
  
  while(it != subSystems.end())
  {
    m_statusModel->addSubSystem(*it, nbProcs);
    it++;
  }
  
  m_networkComms[index]->sendActivateSimulation(nbProcs, hosts.join(","));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::simulationStatus(const QString & subSysName, int rank, 
                                    const QString & status)
{
  m_statusModel->setWorkerStatus(subSysName, rank, status);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::sshError()
{
  QProcess * process = static_cast<QProcess *>(sender());
  QModelIndex index = m_launchServerProcs.key(process);
  
  if(process != NULL && index.isValid())
  {
    m_timers[index]->stop();
    
    QString errorMsg = m_launchServerProcs[index]->readAllStandardError();
    GlobalLog::error(errorMsg);
    
    // stop the process (send SIGKILL-like signal)
    m_launchServerProcs[index]->kill();
  }
  else
  {
    QString errorMsg = "Unexpected call of ClientKernel::sshError() from a "
    "sender of type: %1"; 
    if(sender() != NULL)
      GlobalLog::error(errorMsg.arg(sender()->metaObject()->className()));
    else
      GlobalLog::error(errorMsg.arg("<unknown type>"));
  }
}

/****************************************************************************
 
 PUBLIC SLOTS
 
 ****************************************************************************/

void ClientKernel::connectSimulation(const QModelIndex & index, 
                                     const TSshInformation & info)
{
  if(m_networkComms.contains(index) && m_networkComms.value(index) != NULL 
     && m_networkComms[index]->isConnected())
    GlobalLog::error("This simulation is already connected.");
  else
  {
    ClientNetworkComm * comm = new ClientNetworkComm();
    QTimer * timer = new QTimer(this);
    QProcess * launchProcess = new QProcess(this);
    
    // forward some signals from the network layer to the upper level
    connect(comm, SIGNAL(dirContents(const QString &, const QStringList &, const QStringList &)), 
            this, SIGNAL(dirContents(const QString &, const QStringList &, const QStringList &)));
    
    // connectSig is a macro defined at the top of this file
    connectSig(comm, newTree(const QDomDocument &));
    connectSig(comm, disconnectFromServer());
    connectSig(comm, connected());
    connectSig(comm, ack(COOLFluiD::network::NetworkFrameType));
    connectSig(comm, nack(COOLFluiD::network::NetworkFrameType));
    connectSig(comm, m_abstractTypes(const QStringList &));
    connectSig(comm, concreteTypes(const QStringList &));
    connectSig(comm, ack(COOLFluiD::network::NetworkFrameType));
    connectSig(comm, nack(COOLFluiD::network::NetworkFrameType)); 
    connectSig(comm, m_hostList(const QList<COOLFluiD::network::HostInfos> &));
    connectSig(comm, simulationStatus(const QString &, int, const QString &));
    connectSig(comm, subsystemList(const QStringList &));
    
    if(m_networkComms.contains(index))
      delete m_networkComms.value(index);
    
    if(m_timers.contains(index))
      delete m_timers.value(index);
    
    if(m_launchServerProcs.contains(index))
      delete m_launchServerProcs.value(index);
    
    m_networkComms[index] = comm;
    
    m_timers[index] = timer;
    m_launchServerProcs[index] = launchProcess;
    m_commsSshInfo[index] = info;
    
    connect(timer, SIGNAL(timeout()), this, SLOT(tryToConnect()));
    connect(launchProcess, SIGNAL(readyReadStandardError()), this, SLOT(sshError()));
    
    this->connectToServer(index);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::disconnectSimulation(const QModelIndex & index, bool shutServer)
{
  if(!m_networkComms.contains(index) || m_networkComms.value(index) == NULL 
     || !m_networkComms[index]->isConnected())
    GlobalLog::error("This simulation is not connected.");
  else
    m_networkComms[index]->disconnectFromServer(shutServer);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::runSimulation(const QModelIndex & index)
{
  m_networkComms[index]->sendRunSimulation();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::stopSimulation(const QModelIndex & index)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::activateSimulation(const QModelIndex & index)
{
  if(!m_networkComms.contains(index) || m_networkComms.value(index) == NULL 
     || !m_networkComms[index]->isConnected())
    GlobalLog::error("This simulation is not connected.");
  else
    m_networkComms[index]->sendGetSubSystemList();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::deactivateSimulation(const QModelIndex & index)
{
  if(!m_networkComms.contains(index) || m_networkComms.value(index) == NULL 
     || !m_networkComms[index]->isConnected())
    GlobalLog::error("This simulation is not connected.");
  else
    m_networkComms[index]->sendDeactivateSimulation();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::disconnectFromServer()
{
  QModelIndex index = this->getSimIndex(sender());
  
  if(index.isValid())
  {
    GlobalLog::message("Disconnected from the server.");
    m_treeModel->setSimConnectedStatus(index, false);
  }
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::ack(NetworkFrameType type)
{
  QModelIndex index = this->getSimIndex(sender());
  
  switch(type)
  {
    case NETWORK_OPEN_FILE:
      m_networkComms[index]->sendGetTree();
      break;
      
    case NETWORK_ACTIVATE_SIMULATION:
      m_treeModel->setSimActiveState(index, true);
      break;
      
    case NETWORK_DEACTIVATE_SIMULATION:
      m_treeModel->setSimActiveState(index, false);
      m_statusModel->clear();
      break;
      
    case NETWORK_SIMULATION_RUNNING:
      m_treeModel->setSimReadOnly(index, true);
      break;
      
    case NETWORK_RUN_SIMULATION:
      m_treeModel->setSimReadOnly(index, false);
      break;
      
    default:
      emit acked(type);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::nack(NetworkFrameType type)
{
  switch(type)
  {
    case NETWORK_ACTIVATE_SIMULATION:
      m_statusModel->clear();
      break;
      
    default:
      GlobalLog::message("Unexpected NACK received");
      break;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::addNode(const QString & abstractType)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::renameNode(const QDomNode & node, const QString & newName)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::deleteNode(const QDomNode & node)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientKernel::commitChanges(const QDomDocument & doc)
{
  ClientNetworkComm * comm = m_networkComms[m_treeModel->getCurrentSimulation()];
  comm->sendModifyNode(doc);
}
