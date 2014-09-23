#include <QtCore>

#include "Config/BuilderParser.hh"
#include "Common/StringOps.hh"

#include "ClientServer/server/MPIListener.hh"

using namespace MPI;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::server;

MPIListener::MPIListener(unsigned int waitingTime)
: QThread(),
m_waitingTime(waitingTime)
{
  m_listening = false;
  m_receivingAcksComm = COMM_NULL;
  
  m_receivingAcksMutex = new QMutex();
  m_timer = new QTimer(this);
  
  connect(m_timer, SIGNAL(timeout()), this, SLOT(check()));
  
  // register classes to Qt. This will allow to throw a signal with objects 
  // of these classes.
  qRegisterMetaType<MPI::Intercomm>("MPI::Intercomm"); 
  qRegisterMetaType<COOLFluiD::Config::BuilderParserFrameInfo>("COOLFluiD::Config::BuilderParserFrameInfo"); 
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MPIListener::addComm(const Intercomm & comm)
{
  m_comms[comm] = MPIListeningInfo();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MPIListener::stopListening()
{
  m_listening = false;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MPIListener::run()
{
  QMap<Intercomm, MPIListeningInfo>::iterator it;
  BuilderParserFrameInfo frameInfo;
  
  if(!m_comms.isEmpty())
  {
    m_listening = true;
    
    this->init(); // initialize the listening proess
    m_timer->start(m_waitingTime);
    this->exec(); // run the thread event loop
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MPIListener::WaitingAckResult MPIListener::waitForAcks(const Intercomm & comm, 
                                                       const ManagerWorkerFrameType & ackType,
                                                       unsigned int timeout)
{
  WaitingAckResult result;
  
  if(m_receivingAcksComm == COMM_NULL && comm != COMM_NULL 
     && m_comms.contains(comm))
  {
    int remoteSize = comm.Get_remote_size();
    
    m_ackType = ackType;
    
    m_receivingAcksComm = comm;
    result = SUCCESS;
    
    // fill the vector
    for(int i = 0 ; i < remoteSize ; i++)
      m_receivedAcks.push_back(false);
  }
  else if(comm == COMM_NULL)
    result = FAILURE_NULL_COMM;
  else if(!m_comms.contains(comm))
    result = FAILURE_UNKNOWN_COMM;
  else
    result = FAILURE_ALREADY_WAITING;
  
  return result;
  
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MPIListener::init()
{
  QMap<Intercomm, MPIListeningInfo>::iterator it = m_comms.begin();
  
  // non-blocking receive on all communicators
  while(it != m_comms.end() && m_listening)
  {
    MPIListeningInfo & info = it.value();
    
    if(!info.m_finished && info.m_ready)
    {
      info.m_request = it.key().Irecv(info.m_data, MPIListeningInfo::BUFFER_SIZE, 
                                      CHAR, MPI_ANY_SOURCE, 0);
      info.m_ready = false;
    }
    
    it++;
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MPIListener::check()
{
  QMap<Intercomm, MPIListeningInfo>::iterator it = m_comms.begin();
  BuilderParserFrameInfo frameInfo;
  
  // check all communicators for some new data
  while(it != m_comms.end() && m_listening)
  {
    MPIListeningInfo & info = it.value();
    
    // if data arrived
    if(!info.m_finished && !info.m_ready && info.m_request.Test())
    {
      // parse frame
      if(!BuilderParser::parseFrame(info.m_data, m_rules, frameInfo))
      {
        std::cerr << BuilderParser::getErrorString() << std::endl;
        m_listening = false;
      }
      else
      {
        QMutexLocker lock(m_receivingAcksMutex);
        
        // if we are waiting for acks for ackType frame from this comm 
        if(/*info.comm == this->receivingAcksComm &&*/ frameInfo.frameType == MGR_WKR_ACK 
           && frameInfo.frameAttributes["type"] == m_rules.getTypeName(m_ackType))
        {
          int wkrRank = StringOps::from_str<int>(frameInfo.frameAttributes["workerRank"]);
          QVector<bool>::iterator it = m_receivedAcks.begin();
          
          m_receivedAcks[wkrRank] = true;
          
          bool waitingFinished = true;
          
          while(it != this-> m_receivedAcks.end() && waitingFinished)
            waitingFinished = *(it++);
          
          if(waitingFinished)
          {
            m_receivingAcksComm = COMM_NULL;
            m_ackType = MGR_WKR_NO_TYPE;
            emit ackArrived();
          }
        }
        else
          emit newFrame(it.key(), frameInfo);
        
        it.value().m_ready = true; // ready to do another non-blocking receive      
      }
    }
    
    it++;
  }
  
  if(!m_listening)
  {
    // reinitialize info structures
    it = m_comms.begin();
    while(it != m_comms.end() && m_listening)
    {
      it.value() = MPIListeningInfo();
      it++;
    }
    m_timer->stop();
    this->exit(0);
  }
  else
    this->init();
}
