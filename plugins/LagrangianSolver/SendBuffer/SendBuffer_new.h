#ifndef COOLFluiD_SendBuffer_hh
#define COOLFluiD_SendBuffer_hh

//////////////////////////////////////////////////////////////////////////////
#include "Common/MPI/MPIError.hh"
#include <vector>
#include "Common/COOLFluiD.hh"
#include "Common/PE.hh"
#include "LagrangianSolver/LagrangianSolverModule.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

namespace LagrangianSolver{

template<typename T>
class SendBuffer
{
private:
    T* m_sendBuffer;
    CFint* m_sendCounts;
    CFint* m_sendRanks;
    T* m_sendBufferOrdered;
    CFuint              m_nbProcesses;
    MPI_Comm            m_comm;
    MPI_Datatype        m_MPIdatatype;
    CFuint m_idx;
    CFuint m_bufferSize;

public:
    SendBuffer();
    ~SendBuffer();
    void setupSize(CFuint size);
    bool sincronize(std::vector<T> &recvBuffer, bool isLastPhoton);
    void push_back(const T &a, const CFuint &rank);
    MPI_Datatype getMPIdatatype() const{ return m_MPIdatatype; }
    void setMPIdatatype(const MPI_Datatype MPIdatatype){ m_MPIdatatype = MPIdatatype; }
};

template<typename T>
SendBuffer<T>::SendBuffer():
{
    m_idx= 0;
    m_comm = Common::PE::GetPE().GetCommunicator();
    m_nbProcesses= Common::PE::GetPE().GetProcessorCount();
}

template<typename T>
SendBuffer<T>::~SendBuffer():
{
    delete[] m_sendCounts;
    delete[] m_sendBuffer;
    delete[] m_sendRanks;
    delete[] m_sendBufferOrdered;
}




template<typename T>
void SendBuffer<T>::setupSize(CFuint size){
    m_bufferSize = size;
    m_sendCounts = new CFuint[m_nbProcesses];
    m_sendBuffer = new T[size];
    m_sendRanks = new CFuint[size];
    m_sendBufferOrdered = new T[size]; //is it necessary?
    m_idx=0;
}

template<typename T>
void SendBuffer<T>::push_back( const T &a, const CFuint &rank){
    cf_assert(rank<m_nbProcesses);
    cf_assrt(m_idx< m_bufferSize);
    m_sendBuffer[m_idx]= a;
    m_sendRanks[m_idx] = rank;
    ++m_sendCounts[rank];
}

template<typename T>
bool SendBuffer<T>::sincronize( std::vector<T> &recvBuffer, bool isLastPhoton){

    //get the number of photons to send and the displacements

    CFuint nbPhotonsSend=0;
    std::vector<CFint> displacements(m_nbProcesses);
    for(CFuint i=0; i< m_nbProcesses ; ++i ){
      displacements[i] = nbPhotonsSend;
      nbPhotonsSend += m_sendCounts[i];
    }

    m_sendBufferOrdered.resize(nbPhotonsSend);

    //get the number of photons to receive
    std::vector<CFint> recvCounts(m_nbProcesses), recvDisps(m_nbProcesses);
    MPI_Alltoall(&m_sendCounts[0], 1 , MPI_UNSIGNED, &recvCounts[0],  1 , MPI_UNSIGNED, m_comm);

    CFuint nbPhotonsRecv=0;
    for(CFuint i=0; i< m_nbProcesses ; ++i ){
      recvDisps[i]   = nbPhotonsRecv;
      nbPhotonsRecv += recvCounts[i];
    }

    recvBuffer.resize(nbPhotonsRecv);

    //copy and organize the data into the new buffer
    //TODO: let's look for a way to do it without an extra buffer!
    std::vector<CFint> tempDisps= displacements;
    CFint *tempDisp;
    for( CFuint i=0; i < m_sendBuffer.size(); ++i ){
      tempDisp= &(tempDisps[ m_sendRanks[i] ]);
      m_sendBufferOrdered[ *tempDisp ] = m_sendBuffer[i];
      ++ *tempDisp;
    }

    MPI_Alltoallv(&m_sendBufferOrdered[0], &m_sendCounts[0], &displacements[0],
                  m_MPIdatatype, &recvBuffer[0], &recvCounts[0],
                  &recvDisps[0], m_MPIdatatype, m_comm );

    //clear the sendbuffers
    m_sendBuffer.clear();
    m_sendRanks.clear();
    for(CFuint i=0; i<m_sendCounts.size(); ++i){
        m_sendCounts[i]=0;
    }

    //check finish condition (all buffers have zero size and all partitions have generated all photons)
    CFuint totalPhotonsRecv;
    nbPhotonsRecv += (isLastPhoton)? 0 : 1 ;

    MPI_Allreduce(&nbPhotonsRecv, &totalPhotonsRecv, 1, MPI_UNSIGNED, MPI_SUM, m_comm);
    return (totalPhotonsRecv == 0);
  }  

}

}


#endif
