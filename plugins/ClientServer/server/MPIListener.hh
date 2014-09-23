#ifndef COOLFluiD_server_MPIListener_h
#define COOLFluiD_server_MPIListener_h

////////////////////////////////////////////////////////////////////////////

#include <QMap>
#include <QMutex>
#include <QVector>
#include <QThread>
#include <mpi.h>

#include "Config/BuilderParserFrameInfo.hh"
#include "Config/ManagerWorkerProtocol.hh"
#include "Config/ManagerWorkerFrameType.hh"
#include "ClientServer/server/MPIListeningInfo.hh"

class QTimer;

namespace COOLFluiD
{
  namespace server
  {
    
    ////////////////////////////////////////////////////////////////////////////
    
    /// @brief Listener for MPI messages.
    
    /// @c MPIListener is a thread that continuously listens to one or more 
    /// MPI intercommunicators. If there are more than one intercommunicator,
    /// all communicators are listened simultanuously.
    
    /// @author Quentin Gasper
    
    class MPIListener : public QThread
    {
      Q_OBJECT 
      
    public:
      
      enum WaitingAckResult
      {
        SUCCESS,
        
        FAILURE_ON_NACK,
        
        FAILURE_ON_TIMEOUT,
        
        FAILURE_ALREADY_WAITING,
        
        FAILURE_UNKNOWN_COMM,
        
        FAILURE_NULL_COMM
      };
      
      /// @brief Constructor.
      MPIListener(unsigned int waitingTime = 100);
      
      /// @brief Adds a communicator to listen to.
      
      /// This method can be called during the listening.
      /// @param comm Communicator to add.
      void addComm(const MPI::Intercomm & comm);
      
      /// @brief Stops the listening.
      
      /// Calling this method will exit the thread execution if it is running.
      void stopListening();
      
      /// @brief Reimplements QThread::run() method.
      
      /// If there is at least one communicator to listen to, the process has three
      /// main steps:
      /// @li call IRecv (non-blocking receive) on ready communicator (new ones 
      /// and those that just received data)
      /// @li wait @c #waitingTime ms
      /// @li call check all communicators for new data
      /// When new data arrived, @c #newFrame signal is emitted. @n
      /// This process is repeated as long as @c #stopListening is not called. 
      /// If a new communicator is added, it will be taken in account on the next
      /// iteration.
      void run();
      
      WaitingAckResult waitForAcks(const MPI::Intercomm & comm, 
                                   const COOLFluiD::Config::ManagerWorkerFrameType & ackType,
                                   unsigned int timeout = 1000);
      
      private slots:
      
      void check();
      
    signals:
      
      /// @brief Signal emitted when new data have been received.
      
      /// @param senderComm Communicator from which the data came.
      /// @param frameInfo Parse frame information
      void newFrame(const MPI::Intercomm & senderComm, 
                    const COOLFluiD::Config::BuilderParserFrameInfo & frameInfo);
      
      void ackArrived();
      
    private:
      
      /// @brief Protocol rules
      COOLFluiD::Config::ManagerWorkerProtocol m_rules;
      
      /// @brief Communicators
      
      /// The key is the communicator. The value holds information relative to
      /// the communicator.
      QMap<MPI::Intercomm, MPIListeningInfo> m_comms;
      
      /// @brief Indicates whether the thread is listening.
      
      /// If @c true, the thread is listening, otherwise it is not.
      bool m_listening;
      
      unsigned int m_waitingTime;
      
      QVector<bool> m_receivedAcks;
      
      MPI::Intercomm m_receivingAcksComm;
      
      QMutex * m_receivingAcksMutex;
      
      COOLFluiD::Config::ManagerWorkerFrameType m_ackType;
      
      QTimer * m_timer;
      
      void init();
      
      
    }; // class MPIReceiver
    
    ////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_MPIReceiver_h