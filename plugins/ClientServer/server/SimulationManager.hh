#ifndef COOLFluiD_server_SimulationManager_h
#define COOLFluiD_server_SimulationManager_h

//////////////////////////////////////////////////////////////////////////////

#include <QThread>
#include <QHash>

#include <mpi.h>

#include "Config/ManagerWorkerProtocol.hh"
#include "Config/ManagerWorkerFrameType.hh"
#include "Config/BuilderParserFrameInfo.hh"

#include "ClientServer/server/MPIListener.hh"


//////////////////////////////////////////////////////////////////////////////

class QString;

namespace COOLFluiD
{
  
  namespace server
  {
    
    class MPIListener;
    
    //////////////////////////////////////////////////////////////////////////////
    
    /// @brief Simulation manager
    
    /// This class is a thread.
    
    /// @author Quentin Gasper
    
    class SimulationManager : public QThread
    {
      Q_OBJECT
      
    public:
      
      SimulationManager();
      
      /// @brief Thread exacution code.
      
      /// Reimplements @c QThread::run() method.
      void run();
      
      /// @brief Spawn workers
      
      /// @param subSysName SubSystem name
      /// @param subSysTyoe SubSystem type
      /// @param nbProcs Number of workers to spawn
      /// @param hosts Hosts on which the worker will be spawned
      void spawn(const QString & subSysName, const QString & subSysType, 
                 unsigned int nbProcs, const QString & hosts);
      
      /// @brief Sends a configuration tree to all workers
      
      /// @param configTree Configuration tree to send
      void configure(const XMLNode & configTree);
      
      /// @brief Interconnect workers of different subsystems
      
      /// If one of the provided susystems is either not spawned or does not 
      /// exist, nothing is done. The left subsystem opens a port to allow the
      /// right subsystem to connect.
      /// @param leftRole Left subsystem
      /// @param rightRole Right subsystem
      void connectSubSystems(const QString & leftRole, const QString & rightRole);
      
      /// @brief Sends to all workers a request to exit
      void exitWorkers();
      
      private slots:
      
      void processFrame(MPI::Intercomm m_comm, const COOLFluiD::Config::BuilderParserFrameInfo & frameInfo);
      
      void ackArrived();
      
    signals:
      
      /// @brief Signal emitted when a message comes from a worker
      
      /// @param string The message
      void message(const QString & m_string);
      
      void simulationStatus(const QString & subSysName, int rank, 
                            const QString & status);
      
      void simulationTree(const XMLNode & tree);
      
      void simulationFinished();
      
      void spawned();
      
    private:
      
      /// @brief Subsystem intercommunicators
      
      /// The key is the subsystem name. The value is the associated 
      /// intercommunicator.
      QHash<QString, MPI::Intercomm> m_subSysComms;
      
      /// @brief Manager-Worker m_protocol rules
      COOLFluiD::Config::ManagerWorkerProtocol m_protocol;
      
      MPIListener m_listener;
      
      bool m_spawning;
      
      /// @brief Builds and send a frame to another subsystem
      
      /// @param subSys Subsystem name
      /// @param frameInfo Frame information
      void send(const QString & subSys,
                const COOLFluiD::Config::BuilderParserFrameInfo & frameInfo);
      
      /// @brief Builds and send a frame to the manager.
      
      /// @param frameInfo Frame information
      void send(const COOLFluiD::Config::BuilderParserFrameInfo & frameInfo);
      
      /// @brief Processes an ack
      
      /// @param frameInfo ACK frame information
      void processAck(const COOLFluiD::Config::BuilderParserFrameInfo & frameInfo);
      
      void waitForAcks(MPI::Intercomm comm, COOLFluiD::Config::ManagerWorkerFrameType type);
      
      
    }; // class SimulationManager
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_SimulationManager_h
