#include <QtCore>
#include <QtNetwork>
#include <QtXml>
#include <string>

#include "Common/xmlParser.h"
#include "Common/StringOps.hh"

#include "Config/BuilderParser.hh"
#include "Config/BuilderParserFrameInfo.hh"

#include "ClientServer/client/GlobalLog.hh"

#include "ClientServer/client/ClientNetworkComm.hh"

using namespace COOLFluiD::client;
using namespace COOLFluiD::network;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;

ClientNetworkComm::ClientNetworkComm()
{
  m_socket = new QTcpSocket(this);
  
  // connect useful signals to slots
  connect(m_socket, SIGNAL(readyRead()), this, SLOT(newData()));
  connect(m_socket, SIGNAL(disconnected()), this, SLOT(disconnected()));
  connect(m_socket, SIGNAL(error(QAbstractSocket::SocketError)), this,
          SLOT(socketError(QAbstractSocket::SocketError)));
  
  connect(m_socket, SIGNAL(connected()), this, SLOT(connectionEstablished()));
  
  m_blockSize = 0;
  m_requestDisc = false;
  m_connectedToServer = false;
  m_skipRefused = false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClientNetworkComm::~ClientNetworkComm()
{
  delete m_socket;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendCloseFile() 
{
  if(!this->checkConnected())
    return false;
  
  return this->buildAndSend(NETWORK_CLOSE_FILE);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendGetAbstractTypes(const QString & typeName)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_GET_ABSTRACT_TYPES);
  fi.frameAttributes["typeName"] = typeName.toStdString();
  return this->buildAndSend(fi);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendGetConcreteTypes(const QString & typeName)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_GET_CONCRETE_TYPES);
  fi.frameAttributes["typeName"] = typeName.toStdString();
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendAddNode(const QDomNode & node,
                                    const QString & type,
                                    const QString & absType)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false ;
  
  fi.setFrameType(NETWORK_ADD_NODE);
  fi.frameAttributes["path"] = this->getNodePath(node).toStdString();
  fi.frameAttributes["type"] = type.toStdString();
  fi.frameAttributes["absType"] = absType.toStdString();
  
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendRenameNode(const QDomNode & node,
                                       const QString & newName)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_RENAME_NODE);
  fi.frameAttributes["path"] = this->getNodePath(node).toStdString();
  fi.frameAttributes["newName"] = newName.toStdString();
  
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendDeleteNode(const QDomNode & node)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_DELETE_NODE);
  fi.frameAttributes["path"] = this->getNodePath(node).toStdString();
  
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendModifyNode(const QDomDocument & data)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_MODIFY_NODE);
  fi.frameData = XMLNode::parseString(data.toString().toStdString().c_str());
  
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientNetworkComm::connectToServer(const QString & hostAddress, quint16 port, 
                                        bool skipRefused)
{
  m_skipRefused = skipRefused;
  m_socket->connectToHost(hostAddress, port);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientNetworkComm::disconnectFromServer(bool shutServer)
{
  if(shutServer)
    this->buildAndSend(NETWORK_SHUTDOWN_SERVER);
  
  m_requestDisc = true;
  m_connectedToServer = false;
  
  // close the m_socket
  m_socket->abort();
  m_socket->close();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendOpenFile(const QString & filename)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_OPEN_FILE);
  fi.frameAttributes["filename"] = filename.toStdString();
  
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendOpenDir(const QString & dirname, 
                                    bool includeFiles,
                                    const QStringList & extensions,
                                    bool includeNoExtension)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_READ_DIR);
  
  fi.frameAttributes["dirPath"] = dirname.toStdString();
  fi.frameAttributes["includeFiles"] = includeFiles ? "yes" : "no";
  fi.frameAttributes["extensions"] = extensions.join("*").toStdString();
  fi.frameAttributes["includeNoExtension"] = includeNoExtension ? "yes" : "no";
  
  return this->buildAndSend(fi);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendGetTree()
{
  if(!this->checkConnected())
    return false;
  
  return this->buildAndSend(NETWORK_GET_TREE);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendGetHostList()
{
  if(!this->checkConnected())
    return false;
  
  return this->buildAndSend(NETWORK_GET_HOST_LIST);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::isConnected() const
{
  return m_connectedToServer;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendConfig(const QDomDocument & config) 
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_CONFIGURE);
  fi.frameData = XMLNode::parseString(config.toString().toStdString().c_str());
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendCreateDir(const QString & path, const QString & name)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_CREATE_DIR);
  fi.frameAttributes["path"] = path.toStdString();
  fi.frameAttributes["dirName"] = name.toStdString();
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendSaveConfig(const QString & path, const QDomDocument & config)
{
  BuilderParserFrameInfo fi;
  QDomDocument toSend = config;
  
  if(!this->checkConnected())
    return false ;
  
  if(toSend.firstChild().nodeType() == QDomNode::ProcessingInstructionNode)
    toSend.replaceChild(toSend.childNodes().at(1), toSend.firstChild()); 
  
  fi.setFrameType(NETWORK_SAVE_CONFIG);
  fi.frameAttributes["filename"] = path.toStdString();
  fi.frameData = XMLNode::parseString(toSend.toString().toStdString().c_str());
  return this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendRunSimulation()
{
  if(!this->checkConnected())
    return false;
  
  return this->buildAndSend(NETWORK_RUN_SIMULATION);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendActivateSimulation(int nbProcs, const QString & hosts)
{
  BuilderParserFrameInfo fi;
  
  if(!this->checkConnected())
    return false;
  
  fi.setFrameType(NETWORK_ACTIVATE_SIMULATION);
  fi.frameAttributes["nbProcs"] = QString::number(nbProcs).toStdString();
  fi.frameAttributes["hosts"] = hosts.toStdString();
  
  return this->buildAndSend(fi); 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendDeactivateSimulation()
{
  if(!this->checkConnected())
    return false;
  
  return this->buildAndSend(NETWORK_DEACTIVATE_SIMULATION);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::sendGetSubSystemList()
{
  if(!this->checkConnected())
    return false;
  
  return this->buildAndSend(NETWORK_GET_SUBSYSTEM_LIST);
}

/****************************************************************************
 
 PRIVATE METHOD
 
 ****************************************************************************/

int ClientNetworkComm::send(const QString & frame) const
{
  int charsWritten;
  
  QByteArray block;
  QDataStream out(&block, QIODevice::WriteOnly);
  
  out.setVersion(QDataStream::Qt_4_5); // QDataStream version
  out << (quint32)0;    // reserve 32 bits for the frame data size
  out << frame;
  out.device()->seek(0);  // go back to the beginning of the frame
  out << (quint32)(block.size() - sizeof(quint32)); // write the frame data size
  
  charsWritten = m_socket->write(block);
  m_socket->flush();
  
  return charsWritten;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::buildAndSend(const BuilderParserFrameInfo & frameInfos)
{
  bool success = false;
  std::string frame;
  
  try
  {
    if(!BuilderParser::buildFrame(frameInfos, m_protocol, frame))
      GlobalLog::error(BuilderParser::getErrorString().c_str(), false);
    else
      success = this->send(frame.c_str()) != 0;
  }
  catch(std::string str)
  {
    GlobalLog::error(str.c_str(), false);
  }
  return success;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::buildAndSend(NetworkFrameType type)
{
  BuilderParserFrameInfo fi;
  fi.setFrameType(type);
  this->buildAndSend(fi);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString ClientNetworkComm::getNodePath(const QDomNode & node) const
{
  QDomNode parentNode = node.parentNode();
  QString path;
  
  if(parentNode.isNull()) // if the node has no parent
    return QString();
  else
  {
    path = this->getNodePath(parentNode);
    return path + '/' + node.nodeName();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ClientNetworkComm::checkConnected()
{
  if(!m_connectedToServer)
    GlobalLog::error("Not connected to the server.", false);
  
  return m_connectedToServer;
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void ClientNetworkComm::newData()
{
  BuilderParserFrameInfo fi;
  QString frame;
  QDataStream in(m_socket);
  
  in.setVersion(QDataStream::Qt_4_5); // QDataStream version
  
  // if the server sends two messages very close in time, it is possible that
  // the client never gets the second one.
  // So, it is useful to explicitly read the m_socket until the end is reached.
  while(!m_socket->atEnd())
  {
    // if the data size is not known
    if (m_blockSize == 0)
    {
      // if there are at least 4 bytes to read...
      if (m_socket->bytesAvailable() < (int)sizeof(quint32))
        return;
      
      // ...we read them
      in >> m_blockSize;
    }
    
    if (m_socket->bytesAvailable() < m_blockSize)
      return;
    
    in >> frame;
    
    if(!BuilderParser::parseFrame(frame.toStdString(), m_protocol, fi))
      GlobalLog::error(BuilderParser::getErrorString().c_str(), false);
    else
    {
      switch(fi.frameType)
      {
        case NETWORK_MESSAGE:
          GlobalLog::message(fi.frameAttributes["value"].c_str(), true);
          break;
          
        case NETWORK_ERROR:
          GlobalLog::error(fi.frameAttributes["value"].c_str(), true);
          break;
          
          // if the server sends directory contents
        case NETWORK_DIR_CONTENTS:
        {
          QString strDirs = fi.frameAttributes["dirs"].c_str();
          QString strFiles = fi.frameAttributes["files"].c_str();
          QStringList dirs = strDirs.split("*", QString::SkipEmptyParts);
          QStringList files = strFiles.split("*", QString::SkipEmptyParts);
          
          emit dirContents(fi.frameAttributes["dirPath"].c_str(), dirs, files);
          break;
        }
          
          // if the server sends the tree
        case NETWORK_TREE :
        {
          QDomDocument doc;
          doc.setContent((QString) fi.frameData.createXMLString());
          emit newTree(doc);
          break;
        }
          
          // if the server sends an ACK
        case NETWORK_ACK:
          emit ack(m_protocol.convertToType(fi.frameAttributes["type"].c_str()));
          break;
          
          // if the server sends an NACK
        case NETWORK_NACK:
          emit nack(m_protocol.convertToType(fi.frameAttributes["type"].c_str()));
          break;
          
          // if the server sends the abstract types list
        case NETWORK_ABSTRACT_TYPES :
        {
          QString types = fi.frameAttributes["typeList"].c_str();
          emit abstractTypes(types.split(", "));
          break;
        }
          
          // if the server sends the concrete types list
        case NETWORK_CONCRETE_TYPES :
        {
          QString types = fi.frameAttributes["typeList"].c_str();
          emit concreteTypes(types.split(", "));
          break;
        }
          
          // if the server sends the host list
        case NETWORK_HOST_LIST :
        {
          QList<HostInfos> list;
          QDomDocument doc;
          QDomNodeList childNodes;
          QString err;
          XMLNode & data = fi.frameData;
          
          childNodes = doc.childNodes();
          
          for(int i = 0 ; i < data.nChildNode() ; i++)
          {
            HostInfos hi;
            XMLNode child = data.getChildNode(i);
            
            hi.m_hostname = child.getAttribute("name");
            
            hi.m_nbSlots = QVariant(child.getAttribute("nbSlots")).toInt();
            hi.m_maxSlots = QVariant(child.getAttribute("maxSlots")).toInt();
            
            list.append(hi);
          }
          
          emit hostList(list);
        }
          break;
          
        case NETWORK_SIMULATION_STATUS:
        {
          QString subSysName = fi.frameAttributes["subSysName"].c_str();
          QString status = fi.frameAttributes["value"].c_str();
          int rank = StringOps::from_str<int>(fi.frameAttributes["workerRank"]);
          
          emit simulationStatus(subSysName, rank, status);
        }
          break;
          
        case NETWORK_SUBSYSTEM_LIST:
          emit subsystemList(QString(fi.frameAttributes["subSystems"].c_str()).split(" "));
          break;
      }
    }
    m_blockSize = 0;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientNetworkComm::connectionEstablished()
{
  m_connectedToServer = true;
  emit connected();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientNetworkComm::disconnected()
{
  if(!m_requestDisc)
  {
    GlobalLog::error("The connection has been closed", false);
  }
  emit disconnectFromServer();
  
  m_connectedToServer = false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ClientNetworkComm::socketError(QAbstractSocket::SocketError err)
{
  if(m_requestDisc)
    return;
  
  if(m_connectedToServer)
  {
    m_requestDisc = true;
    m_socket->disconnectFromHost();
  }
  
  switch (err)
  {
    case QAbstractSocket::RemoteHostClosedError:
      GlobalLog::error("Remote connection closed", false);
      break;
      
    case QAbstractSocket::HostNotFoundError:
      GlobalLog::error("Host was not found", false);
      break;
      
    case QAbstractSocket::ConnectionRefusedError:
      if(!m_skipRefused)
        GlobalLog::error("Connection refused. Please check if the server is "
                         "running.", false);
      else
        m_skipRefused = false;
      break;
      
    default:
      GlobalLog::error(QString("The following error occurred: ") +
                       m_socket->errorString(), false);
  }
}
