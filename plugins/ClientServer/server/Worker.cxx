#include <iostream>
#include <string>

#include <QCoreApplication>
#include <QDebug>
#include <mpi.h>

#include "ClientServer/network/FrameInfos.h"
#include "ClientServer/network/FrameBuilderParser.h"
#include "server/MgrWorkerProtocol.h"
#include "server/WorkerThread.h"
#include "server/Worker.h"

using namespace std;
using namespace MPI;
using namespace COOLFluiD::server;
using namespace COOLFluiD::network;

Worker::Worker()
{
 this->m_rank = COMM_WORLD.Get_rank();
 this->m_parent = COMM_WORLD.Get_parent();
 this->m_number = 0;
 this->m_sum = 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Worker::simulate(MgrWorkerProtocol & rules)
{
 Request request = this->m_parent.Irecv(this->m_string, sizeof(m_string), CHAR, 0, 0);
 bool stop = false;
  
 cout << this->getHeader().toStdString() << "Simulating the simulation" << endl;
 
 this->m_parent.Barrier();
 
 for(int i = 0 ; i < 2000000 && !stop ; i++)
 {
  this->m_number = i;
  int localSum = 0;
  int remoteSum = 0;
  QHash<QString, Intercomm>::iterator it = this->m_global.begin();
  
  COMM_WORLD.Allreduce(&this->m_number, &localSum, 1, INT, SUM);
  this->m_sum = localSum;
  
  while(it != this->m_global.end())
  {
   it.value().Allreduce(&this->m_number, &remoteSum, 1, INT, SUM);
   this->m_sum += remoteSum;
   it++;
  }
  
  // every 10 iterations we check for a string to capitalize
  if(i % 10 == 0)
  {
   if(request.Test()) // if a string has arrived
   {
    FrameInfos fi;
    
    if(!FrameBuilderParser::parseFrame(this->m_string, fi, rules))
    {
     cerr << FrameBuilderParser::getStaticError().toStdString() << endl;
     cerr << "Abording" << endl;
     return;
    }
    
    stop = fi.frameType == MGR_WKR_QUIT;
    
    if(fi.frameType == MGR_WKR_STRING)
    {
     QString strFromListener = this->getString(fi.frameAttributes["value"]);
     QString frame;
     
     fi.setFrameType(MGR_WKR_STRING);
     fi.frameAttributes["value"]  = strFromListener;
     FrameBuilderParser::buildFrame(frame, fi, rules);

     this->m_parent.Send(frame.toStdString().c_str(), frame.length() + 1, CHAR, 0, 0);
     
     m_parent.Barrier(); // wait for parent to call Barrier() on its communicator
     
     // do the next non-blocking receive
     request = this->m_parent.Irecv(this->m_string, sizeof(m_string), CHAR, 0, 0);
    }
   } // for "if(request.Test())"
   
  } // for "if(i % 10 == 0)"
  
 } // for "for(...)"
 
 std::m_string header = this->getHeader().toStdString();
 std::m_string dataString = this->getString().toStdString();
 
 if(stop)
  cout << header << "Stopped! Data: " << dataString << endl;
 else
  cout << header << "Loop finished! Data: " << dataString << endl;
 
 // ending barrier
 QHash<QString, Intercomm>::iterator it = this->m_global.begin();
 while(it != this->m_global.end())
 {
  it.value().Barrier();
  it++;
 }
 
 m_parent.Barrier();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int Worker::getNumber() 
{ 
 return this->m_number;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int Worker::getSum() 
{
 return this->m_sum;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString Worker::getString(const QString & string) 
{
 QString infos = "%1 [%2][%3][%4]";
 
 infos = infos.arg(this->m_role).arg(this->m_rank).arg(this->m_number).arg(this->m_sum);
 
 if(!string.isEmpty())
  infos.prepend(string.toUpper() + " -- ");
 
 return infos;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString Worker::getHeader() const
{
 if(m_role.isEmpty())
  return QString("Worker[%1] ").arg(this->m_rank);
 
 return QString("Worker[%1] (%2) ").arg(this->m_rank).arg(this->m_role);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Worker::openPort(const QString & remoteRole)
{
 char myPort[MPI_MAX_PORT_NAME];
 Intercomm remoteComm;
 
 if(this->m_rank == 0)
 {
  Open_port(MPI_INFO_NULL, myPort); // open a port and build the port name

  // send the port name to the manager
  this->m_parent.Send(myPort, sizeof(myPort), CHAR, 0, 0);
 }
  
 // broadcast the port name to all "buddy workers"
 // if this->rank == 0, the string is sent; otherwise the string is received
 COMM_WORLD.Bcast(myPort, sizeof(myPort), CHAR, 0);
  
 m_parent.Barrier();
  
 remoteComm = COMM_WORLD.Accept(myPort, MPI_INFO_NULL, 0);
 
 cout << this->getHeader().toStdString() << "Found " << remoteComm.Get_remote_size() 
   << " \"" << remoteRole.toStdString() << "\" process(es)." << endl;
 
 remoteComm.Barrier();
 COMM_WORLD.Barrier();
 
 this->m_global[remoteRole] = remoteComm;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Worker::connectToPort(const QString & portName, const QString & remoteRole)
{
 Intercomm remoteComm;
 
 remoteComm = COMM_WORLD.Connect(portName.toStdString().c_str(), MPI_INFO_NULL, 0);
 
 cout << this->getHeader().toStdString() << "Found " << remoteComm.Get_remote_size() 
   << " \"" << remoteRole.toStdString() << "\" process(es)." << endl;
 
 remoteComm.Barrier();
 COMM_WORLD.Barrier();
 
 this->m_global[remoteRole] = remoteComm;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Worker::setRole(const QString & role)
{
 this->role = role;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString Worker::getRole() const
{
 return this->m_role;
}
