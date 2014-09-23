#include <string>

#include <QtCore>

#include "Common/PE.hh"
#include "Common/StringOps.hh"
#include "Common/xmlParser.h"
#include "Config/BuilderParser.hh"
#include "Config/BuilderParserFrameInfo.hh"
#include "Config/ManagerWorkerFrameType.hh"

#include "ClientServer/server/MPIReceiver.hh"
#include "ClientServer/server/SimulationManager.hh"

using namespace MPI;
using namespace std;

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

using namespace COOLFluiD::server;

SimulationManager::SimulationManager()
{
  connect(&m_listener, SIGNAL(newFrame(const MPI::Intercomm &, const COOLFluiD::Config::BuilderParserFrameInfo &)), 
          this, SLOT(processFrame(const MPI::Intercomm &, const COOLFluiD::Config::BuilderParserFrameInfo &)));
  connect(&m_listener, SIGNAL(ackArrived()), this, SLOT(ackArrived()));
  m_spawning = false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::run()
{
  m_listener.start();
  //  while(1)
  //  {
  //   BuilderParserFrameInfo frameInfo;
  //   Intercomm comm = this->subSysComms.begin().value();
  //   
  //   MPIReceiver::receive(frameInfo, comm, m_protocol);
  //   this->processFrame(comm, frameInfo);
  //  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::spawn(const QString & subSysName, const QString & subSysType, 
                              unsigned int nbProcs, const QString & hosts)
{
  BuilderParserFrameInfo frameInfo;
  
  frameInfo.setFrameType(MGR_WKR_SET_SUBSYS);
  frameInfo.frameAttributes["name"] = subSysName.toStdString();
  frameInfo.frameAttributes["type"] = subSysType.toStdString();
  
  if(!subSysName.isEmpty() && !m_subSysComms.contains(subSysName))
  {
    
    for(int i = 0 ; i < nbProcs ; i++)
      emit simulationStatus(subSysName, i, WorkerStatus::Convert::to_str(WorkerStatus::STARTING).c_str()); 
    
    m_subSysComms[subSysName] = PE::spawn(nbProcs, hosts.toStdString().c_str());
    m_subSysComms[subSysName].Barrier();
    m_listener.addComm(m_subSysComms[subSysName]);
    m_spawning = true;
    m_listener.waitForAcks(m_subSysComms[subSysName], MGR_WKR_SET_SUBSYS);
    this->send(subSysName, frameInfo);
    //   this->waitForAcks(this->subSysComms[subSysName], MGR_WKR_SET_SUBSYS);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::configure(const XMLNode & configTree)
{
  BuilderParserFrameInfo frameInfo;
  QHash<QString, Intercomm>::iterator it = m_subSysComms.begin();
  
  frameInfo.setFrameType(MGR_WKR_CONFIGURE);
  frameInfo.frameData = configTree;
  
  // send to all subsystem groups
  while(it != m_subSysComms.end())
  {
    this->send(it.key(), frameInfo);
    it++;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::connectSubSystems(const QString & leftRole, const QString & rightRole)
{
  char data[65536];
  BuilderParserFrameInfo fi;
  string portName;
  Intercomm leftComm;
  Intercomm rightComm;
  
  if(m_subSysComms.contains(leftRole) && m_subSysComms.contains(rightRole))
  {
    leftComm = m_subSysComms[leftRole];
    rightComm = m_subSysComms[rightRole];
    
    emit message("connecting..." + leftRole + " " + rightRole);
    
    // ask the group to open a port
    fi.setFrameType(MGR_WKR_OPEN_PORT);
    fi.frameAttributes["remoteRole"] = rightRole.toStdString();
    this->send(leftRole, fi);
    
    // wait for the group to send a port name
    MPIReceiver::receive(fi, leftComm, m_protocol);
    
    if(fi.frameType == MGR_WKR_PORT_NAME)
    {
      portName = fi.frameAttributes["value"];
      
      std::cout << "Port name: " << portName << std::endl;
      
      // // //    leftComm.Barrier();
      
      // send the port name to the new group and ask it to connect
      fi.setFrameType(MGR_WKR_CONNECT);
      fi.frameAttributes["portName"] = portName;
      fi.frameAttributes["remoteRole"] = leftRole.toStdString();
      
      this->send(rightRole, fi);
      
      // // //    rightComm.Barrier();
      // // //    leftComm.Barrier();
    }
  }
}

/****************************************************************************
 
 PRIVATE METHODS
 
 ****************************************************************************/

void SimulationManager::send(const QString & subSys,
                             const BuilderParserFrameInfo & frameInfo)
{
  string frame;
  
  if(!BuilderParser::buildFrame(frameInfo, m_protocol, frame))
    qDebug() << BuilderParser::getErrorString().c_str();
  else
  {
    Intercomm comm = m_subSysComms[subSys];
    const char * data = frame.c_str();
    
    for(int i = 0 ; i < comm.Get_remote_size() ; i++)
      comm.Send(data, frame.length() + 1, CHAR, i, 0);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::send(const BuilderParserFrameInfo & frameInfo)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::processAck(const BuilderParserFrameInfo & frameInfo)
{
  BuilderParserFrameInfo fi;
  string strRank = frameInfo.frameAttributes.find("workerRank")->second;
  string strType = frameInfo.frameAttributes.find("type")->second;
  int workerRank = QVariant(strRank.c_str()).toInt();
  unsigned int type = m_protocol.getFrameType(strType);
  
  switch(type)
  {
    case MGR_WKR_CONFIGURE:
      fi.setFrameType(MGR_WKR_SIMULATE);
      this->send("SubSystem", fi);
      break;
      
    case MGR_WKR_SIMULATE:
      emit simulationFinished();
      break;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::exitWorkers()
{
  BuilderParserFrameInfo fi;
  fi.setFrameType(MGR_WKR_EXIT);
  
  m_listener.waitForAcks(m_subSysComms["SubSystem"], MGR_WKR_EXIT);
  
  this->send("SubSystem", fi);
  //  this->waitForAcks(this->subSysComms["SubSystem"], MGR_WKR_EXIT);
  m_subSysComms["SubSystem"].Barrier();
  m_subSysComms.remove("SubSystem");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::processFrame(Intercomm comm, 
                                     const BuilderParserFrameInfo & frameInfo)
{
  switch(frameInfo.frameType)
  {
    case MGR_WKR_ACK:
      this->processAck(frameInfo);
      break;
      
    case MGR_WKR_STRING:
      emit message(frameInfo.getAttribute("value").c_str());
      break;
      
    case MGR_WKR_STATUS:
    {
      QString subSysName = frameInfo.getAttribute("subSysName").c_str();
      QString status = frameInfo.getAttribute("value").c_str();
      int rank = StringOps::from_str<int>(frameInfo.getAttribute("workerRank"));
      
      if(subSysName.isEmpty())
        subSysName = m_subSysComms.key(comm);
      
      emit simulationStatus(subSysName, rank, status);
    }
      break;
      
    case MGR_WKR_TREE:
      emit simulationTree(frameInfo.frameData);
      break;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::waitForAcks(Intercomm comm, ManagerWorkerFrameType type)
{
  int nbWorkers = comm.Get_remote_size();
  bool receivedAcks[nbWorkers];
  int receivedAcksCounter = 0;
  BuilderParserFrameInfo frameInfo;
  int wkrRank;
  std::string typeStr = m_protocol.getTypeName(type);
  
  for(int i = 0 ; i < nbWorkers ; i++)
    receivedAcks[i] = false;
  
  while(receivedAcksCounter < nbWorkers)
  {
    MPIReceiver::receive(frameInfo, comm, m_protocol);
    // if this is an ACK for the expected type
    if(frameInfo.frameType == MGR_WKR_ACK && 
       frameInfo.frameAttributes["type"] == typeStr)
    {
      wkrRank = StringOps::from_str<int>(frameInfo.frameAttributes["workerRank"]);
      
      if(!receivedAcks[wkrRank])
      {
        receivedAcks[wkrRank] = true;
        receivedAcksCounter++;
      }
      
    }
    else // we process the frame as a normal one
      this->processFrame(comm, frameInfo);  
  } 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SimulationManager::ackArrived()
{
  if(m_spawning)
  {
    m_spawning = false;
    emit spawned();
  }
}
