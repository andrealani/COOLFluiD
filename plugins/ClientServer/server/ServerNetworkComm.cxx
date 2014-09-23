#include <iostream>
#include <string>

#include <QtNetwork>
#include <QtXml>
#include <QtCore>

#include "Common/xmlParser.h"
#include "Common/StringOps.hh"

#include "Config/ConfigArgs.hh"
#include "Config/ConverterTools.hh"
#include "Config/BuilderParser.hh"
#include "Config/BuilderParserFrameInfo.hh"

#include "ClientServer/network/NetworkException.hh"

#include "ClientServer/server/UnknownClientIdException.hh"
#include "ClientServer/server/ServerSimulation.hh"

#include "ClientServer/server/ServerNetworkComm.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::network;
using namespace COOLFluiD::server;

ServerNetworkComm::ServerNetworkComm(QString hostAddress, quint16 port)
{
  bool local = hostAddress == "127.0.0.1";
  
  m_server = new QTcpServer(this);
  
  if(!local)
    m_localSocket = new QTcpServer(this);
  else
    m_localSocket = NULL;
  
  if(!m_server->listen(QHostAddress(hostAddress), port))
  {
    QString message = QString("Cannot listen %1 on port %2 : %3")
    .arg(hostAddress)
    .arg(port)
    .arg(m_server->errorString());
    throw NetworkException(FromHere(), message.toStdString());
  }
  
  if(!local && !m_localSocket->listen(QHostAddress("127.0.0.1"), port))
  {
    QString message = QString("Cannot listen 127.0.0.1 on port %2 : %3")
    .arg(port)
    .arg(m_server->errorString());
    throw NetworkException(FromHere(), message.toStdString());
  }
  
  m_lastClientId = 0;
  
  connect(m_server, SIGNAL(newConnection()), this, SLOT(newClient()));
  
  if(!local)
    connect(m_localSocket, SIGNAL(newConnection()), this, 
            SLOT(newClient()));
  
  m_server->setMaxPendingConnections (1);
  m_blockSize = 0;
  m_bytesRecieved = 0;
  m_bytesSent = 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ServerNetworkComm::~ServerNetworkComm()
{
  QHash<QTcpSocket*, QDomNode>::iterator it = m_clients.begin();
  
  while(it != m_clients.end())
  {
    delete it.key();
    it++;
  }
  
  m_server->close();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerNetworkComm::buildAndSend(QTcpSocket * client, const BuilderParserFrameInfo & frameInfos)
{
  bool success = false;
  std::string frame;
  
  if(!BuilderParser::buildFrame(frameInfos, m_networkProtocol, frame))
    std::cerr << "Could not build the frame: " << BuilderParser::getErrorString() << frameInfos.frameType << std::endl;
  else
    success = this->send(client, frame.c_str()) != 0;
  
  return success;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int ServerNetworkComm::send(QTcpSocket * client, const QString & frame)
{
  QByteArray block;
  QDataStream out(&block, QIODevice::WriteOnly);
  int count = 0; // total bytes sent
  
  out.setVersion(QDataStream::Qt_4_5);
  // reserving 2 bytes to store the data size 
  // (frame size without these 2 bytes)
  out << (quint32)0;
  out << frame;
  out.device()->seek(0); // go back to the beginning of the frame
  out << (quint32)(block.size() - sizeof(quint32)); // store the data size
  
  if(client == NULL)
  {
    QHash<QTcpSocket *, QDomNode>::iterator it = m_clients.begin();
    
    while(it != m_clients.end())
    {
      client = it.key();
      count += client->write(block);
      m_bytesSent += count;
      client->flush();
      it++;
    }
  }
  
  else
  {
    count = client->write(block);
    m_bytesSent += count;
    
    client->flush();
  }
  
  return count;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerNetworkComm::sendDirContents(int clientId, const QString & path,
                                        const QStringList & dirs, 
                                        const QStringList & files)
{
  bool success = false;
  
  try
  {
    BuilderParserFrameInfo fi;
    QTcpSocket * socket = this->getSocket(clientId);
    
    fi.setFrameType(NETWORK_DIR_CONTENTS);
    fi.frameAttributes["dirPath"] = path.toStdString();
    fi.frameAttributes["dirs"] = dirs.join("*").toStdString();
    fi.frameAttributes["files"] = files.join("*").toStdString();
    
    success = this->buildAndSend(socket, fi);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
  
  return false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerNetworkComm::sendAbstractTypes(int clientId, const QString & typeName,
                                          const QStringList & typeList)
{
  bool success = false;
  
  try
  {
    BuilderParserFrameInfo fi;
    QTcpSocket * socket = this->getSocket(clientId);
    
    fi.setFrameType(NETWORK_ABSTRACT_TYPES);
    fi.frameAttributes["typeName"] = typeName.toStdString();
    fi.frameAttributes["typeList"] = typeList.join(", ").toStdString();
    
    success = this->buildAndSend(socket, fi);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
  
  return false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerNetworkComm::sendConcreteTypes(int clientId, const QString & typeName,
                                          const QStringList & typeList)
{
  bool success = false;
  
  try
  {
    BuilderParserFrameInfo fi;
    QTcpSocket * socket = this->getSocket(clientId);
    
    fi.setFrameType(NETWORK_CONCRETE_TYPES);
    fi.frameAttributes["typeName"] = typeName.toStdString();
    fi.frameAttributes["typeList"] = typeList.join(", ").toStdString();
    
    
    success = this->buildAndSend(socket, fi);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
  
  return false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerNetworkComm::sendHostList(int clientId, 
                                     const QList<COOLFluiD::network::HostInfos> & hosts)
{
  bool success = false;
  QList<COOLFluiD::network::HostInfos>::const_iterator it = hosts.begin();
  
  try
  {
    BuilderParserFrameInfo fi;
    QTcpSocket * socket = this->getSocket(clientId);
    QDomDocument doc;
    
    fi.setFrameType(NETWORK_HOST_LIST);
    
    while(it != hosts.end())
    {
      QDomElement node = doc.createElement("host");
      
      node.setAttribute("name", it->m_hostname);
      node.setAttribute("nbSlots", it->m_nbSlots);
      
      if(it->m_maxSlots != 0)
        node.setAttribute("maxSlots", it->m_maxSlots);
      
      doc.appendChild(node);
      it++;
    }
    
    fi.frameData = XMLNode::parseString(doc.toString().toStdString().c_str());
    
    success = this->buildAndSend(socket, fi);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
  
  return false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendTree(int clientId, const QDomDocument & tree) 
{
  this->sendTree(clientId, XMLNode::parseString(tree.toString().toStdString().c_str()));
  //  try
  //  {
  //   QTcpSocket * m_socket = this->getSocket(clientId);
  //   BuilderParserFrameInfo fi;
  //   
  //   fi.setFrameType(NETWORK_TREE);
  //   fi.frameData = XMLNode::parseString(tree.toString().toStdString().c_str());;
  //   
  //   this->buildAndSend(m_socket, fi);
  //   
  //  }
  //  catch(UnknownClientIdException e)
  //  {
  //   qDebug() << e.what();
  //  }
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendTree(int clientId, const XMLNode & tree) 
{
  try
  {
    QTcpSocket * socket = this->getSocket(clientId);
    BuilderParserFrameInfo fi;
    
    fi.setFrameType(NETWORK_TREE);
    fi.frameData = tree;
    
    this->buildAndSend(socket, fi);
    
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::send(int clientId, const QString & frame) 
{ 
  try
  {
    QTcpSocket * socket = this->getSocket(clientId);
    
    this->send(socket, frame);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendError(int clientId, const QString & message) 
{
  try
  {
    QTcpSocket * socket = this->getSocket(clientId);
    this->sendError(socket, message);
    
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendAck(int clientId, bool success, NetworkFrameType type) 
{
  try
  {
    QTcpSocket * socket = this->getSocket(clientId);
    this->sendAck(socket, success, type);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendMessage(int clientId, const QString & message) 
{
  try
  {
    QTcpSocket * socket = this->getSocket(clientId);
    
    this->sendMessage(socket, message);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendStatus(int clientId,const QString & subSysName, 
                                   int rank, const QString & status)
{
  try
  {
    QTcpSocket * socket = this->getSocket(clientId);
    BuilderParserFrameInfo fi;
    
    fi.setFrameType(NETWORK_SIMULATION_STATUS);
    fi.frameAttributes["subSysName"] = subSysName.toStdString();
    fi.frameAttributes["workerRank"] = StringOps::to_str<int>(rank);
    fi.frameAttributes["value"] = status.toStdString();
    
    this->buildAndSend(socket, fi);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendSubSystemList(int clientId, 
                                          const QStringList & subSystemList)
{
  try
  {
    QTcpSocket * socket = this->getSocket(clientId);
    BuilderParserFrameInfo fi;
    
    fi.setFrameType(NETWORK_SUBSYSTEM_LIST);
    fi.frameAttributes["subSystems"] = subSystemList.join(" ").toStdString();
    
    this->buildAndSend(socket, fi);
  }
  catch(UnknownClientIdException e)
  {
    qDebug() << e.what();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ServerNetworkComm::sendMessage(QTcpSocket * client, const QString & message)
{
  BuilderParserFrameInfo fi;
  
  fi.setFrameType(NETWORK_MESSAGE);
  fi.frameAttributes["value"] = message.toStdString();
  
  return this->buildAndSend(client, fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int ServerNetworkComm::getBytesRecieved() const
{
  return m_bytesRecieved;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int ServerNetworkComm::getBytesSent() const
{
  return m_bytesSent;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendError(QTcpSocket * client, const QString & message)
{
  BuilderParserFrameInfo fi;
  
  fi.setFrameType(NETWORK_ERROR);
  fi.frameAttributes["value"] = message.toStdString();
  
  this->buildAndSend(client, fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::sendAck(QTcpSocket * client, bool success, NetworkFrameType type)
{
  BuilderParserFrameInfo fi;
  
  if(success)
    fi.setFrameType(NETWORK_ACK);
  else
    fi.setFrameType(NETWORK_NACK);
  
  fi.frameAttributes["type"] = m_networkProtocol.getTypeName(type);
  
  this->buildAndSend(client, fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QTcpSocket * ServerNetworkComm::getSocket(int clientId) const
{
  QTcpSocket * socket = NULL;
  
  if(clientId != -1)
  { 
    if(m_clientIds.contains(clientId))
      socket = m_clientIds[clientId];
    
    else
      throw UnknownClientIdException(FromHere(), 
                                     QString("Unknown client id: %1")
                                     .arg(clientId).toStdString());
  }
  
  return socket;
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void ServerNetworkComm::newClient()
{
  QTcpSocket * socket;
  
  socket = m_server->nextPendingConnection();
  
  if(socket == NULL)
    socket = m_localSocket->nextPendingConnection();
  
  // connect useful signals to slots
  connect(socket, SIGNAL(disconnected()), this, SLOT(clientDisconnected()));
  connect(socket, SIGNAL(readyRead()), this, SLOT(newData()));
  
  std::cout << "A new client is connected" << std::endl; 
  
  m_clients[socket] = QDomNode();
  m_clientIds[m_lastClientId] = socket;
  
  emit newClient(m_lastClientId);
  
  m_lastClientId++;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::newData()
{
  // which client has sent data ?
  QTcpSocket * socket = qobject_cast<QTcpSocket *>(sender());
  BuilderParserFrameInfo fi;
  
  int clientId;
  
  QString frame;
  QDataStream in(socket);
  in.setVersion(QDataStream::Qt_4_5);
  
  // if the client sends two messages very close in time, it is possible that
  // the server never gets the second one.
  // So, it is useful to explicitly read the m_socket until the end is reached.
  while(!socket->atEnd())
  {
    // if the data size is not known
    if (m_blockSize == 0)
    {
      // if there are at least 4 bytes to read...
      if (socket->bytesAvailable() < (int)sizeof(quint32))
        return;
      
      // ...we read them
      in >> m_blockSize;
    }
    
    if (socket->bytesAvailable() < m_blockSize)
      return;
    
    in >> frame;
    
    m_bytesRecieved += m_blockSize + (int)sizeof(quint32);
    
    // if parse() returns false, the document is not valid
    if(!BuilderParser::parseFrame(frame.toStdString(), m_networkProtocol, fi))
    {
      std::string error = BuilderParser::getErrorString();
      
      // if error is empty, the document is not a well-formed XML document
      if(error.empty())
        error = "not well-formed document.";
      
      this->sendError(socket, QString("XML parsing error : %1\non string \"%2\"")
                      .arg(error.c_str())
                      .arg(frame));
      
    }
    
    else
    {
      int clientId = m_clientIds.key(socket);
      
      switch(fi.frameType)
      { 
        case NETWORK_READ_DIR:
        {
          QString exts = fi.frameAttributes["extensions"].c_str();
          
          emit openDirectory(clientId, fi.frameAttributes["dirPath"].c_str(), 
                             exts.split("*", QString::SkipEmptyParts),
                             fi.frameAttributes["includeFiles"] == "yes",
                             fi.frameAttributes["includeNoExtension"] == "yes");
          break;
        }
          
          // if the client wants to open a file
        case NETWORK_OPEN_FILE :
          emit openFile(clientId, fi.frameAttributes["filename"].c_str());
          break;
          
          // if the client wants the tree
        case NETWORK_GET_TREE :
          emit getTree(clientId);
          break;
          
          // if the client wants to run the simulation
        case NETWORK_RUN_SIMULATION :
          emit runSimulation(clientId);
          break;
          
          // if the client wants to activate the simulation
        case NETWORK_ACTIVATE_SIMULATION :
        {
          unsigned int nbProcs;
          nbProcs = StringOps::from_str<unsigned int>(fi.frameAttributes["nbProcs"]);
          emit activateSimulation(clientId, nbProcs, fi.frameAttributes["hosts"].c_str());
        }
          break;
          
          // if the client wants to deactivate the simulation
        case NETWORK_DEACTIVATE_SIMULATION :
          emit deactivateSimulation(clientId);
          break;
          
          
          // if the client wants to configure the simulator
        case NETWORK_CONFIGURE :
        {
          QDomDocument doc;
          doc.setContent((QString) fi.frameData.createXMLString());
          emit configureSimulator(clientId, doc);
          break;
        }
          
          // if the client wants to create a new directory
        case NETWORK_CREATE_DIR:
          emit createDirectory(clientId, fi.frameAttributes["path"].c_str(), 
                               fi.frameAttributes["dirName"].c_str());
          break;
          
          // if the client wants to shut the server down
        case NETWORK_SHUTDOWN_SERVER :
          emit shutdownServer(clientId);
          break;
          
          // if the client wants to close the file
        case NETWORK_CLOSE_FILE :
          emit closeFile(clientId);
          break;
          
          // if the client wants the abstract types
        case NETWORK_GET_ABSTRACT_TYPES :
          emit getAbstractTypes(clientId, fi.frameAttributes["typeName"].c_str());
          break;
          
          // if the client wants the concrete types
        case NETWORK_GET_CONCRETE_TYPES :
          emit getConcreteTypes(clientId, fi.frameAttributes["typeName"].c_str());
          break;
          
          // if the client wants to save a configuration on server side
        case NETWORK_SAVE_CONFIG:
        {
          QDomDocument doc;
          doc.setContent((QString) fi.frameData.createXMLString());
          
          emit saveConfiguration(clientId, fi.frameAttributes["filename"].c_str(), doc);
          break;
        }
          
          
          // if the client requests to add a node
        case NETWORK_ADD_NODE :
          emit addNode(clientId, fi.frameAttributes["path"].c_str(), 
                       fi.frameAttributes["nodeName"].c_str(),
                       fi.frameAttributes["type"].c_str(), 
                       fi.frameAttributes["absType"].c_str());
          break;
          
          // if the client requests to rename a node
        case NETWORK_RENAME_NODE :
          emit renameNode(clientId, fi.frameAttributes["path"].c_str(), 
                          fi.frameAttributes["newName"].c_str());
          break;
          
          // if the client requests to delete a node
        case NETWORK_DELETE_NODE :
          emit deleteNode(clientId, fi.frameAttributes["path"].c_str());
          break;
          
          // if the client requests to modify a node
        case NETWORK_MODIFY_NODE :
        {
          QDomDocument doc;
          doc.setContent((QString) fi.frameData.createXMLString());
          emit modifyNode(clientId, doc);
          break;
        }
          
        case NETWORK_GET_HOST_LIST :
          emit getHostList(clientId);
          break;
          
        case NETWORK_GET_SUBSYSTEM_LIST:
          emit getSubSysList(clientId);
          break;
      }
    }
    
    m_blockSize = 0;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::clientDisconnected()
{
  // which client has been disconnected ?
  QTcpSocket * socket = qobject_cast<QTcpSocket *>(sender());
  
  if(socket != NULL)
  {
    int clientId = m_clientIds.key(socket);
    m_clientIds.remove(clientId);
    m_clients.remove(socket);
    
    std::cout << "A client has gone (" << m_clients.size() << " left)\n";
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::message(const QString & message)
{
  this->sendMessage((QTcpSocket*)NULL, message);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ServerNetworkComm::error(const QString & message)
{
  this->sendError((QTcpSocket*)NULL, message);
}
