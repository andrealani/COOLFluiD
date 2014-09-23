#include "ClientServer/server/MPIListeningInfo.hh"

using namespace MPI;
using namespace COOLFluiD::server;

MPIListeningInfo::MPIListeningInfo()
{
  m_counter = 0;
  m_ready = true;
  m_data = new char[BUFFER_SIZE];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MPIListeningInfo::setComm(Intercomm comm)
{
  m_comm = comm;
  m_processCount = comm.Get_remote_size();
  m_finished = m_counter == m_processCount;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MPIListeningInfo::incCounter()
{
  if(!m_finished)
  {
    m_counter++;
    m_finished = m_counter == m_processCount;
  }
}
