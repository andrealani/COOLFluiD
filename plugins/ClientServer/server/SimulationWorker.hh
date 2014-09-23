#ifndef COOLFluiD_server_SimulationWorker_h
#define COOLFluiD_server_SimulationWorker_h

//////////////////////////////////////////////////////////////////////////////

#include <QObject>
#include <QHash>
#include <mpi.h>

#include "Common/PE.hh"

#include "Config/ManagerWorkerProtocol.hh"
#include "Config/BuilderParserFrameInfo.hh"
#include "Config/ManagerWorkerFrameType.hh"

#include "ClientServer/server/MPIListener.hh"
#include "ClientServer/server/ServerSimulation.hh"

namespace COOLFluiD
{
  
  namespace server
  {
    
    class ServerSimulation;
    class MPIListener;
    
    //////////////////////////////////////////////////////////////////////////////
    
    /// @brief Simulation worker
    
    /// @author Quentin Gasper
    
    class SimulationWorker : public QObject
    {
      Q_OBJECT 
      
    public:
      
      /// @brief Constructor
      SimulationWorker();
      
      /// @brief Begins to listen to MPI communicators
      void listen();
      
      QString getLineHeader() const;
      
      public slots:
      
      /// @brief Sends a message to the manager
      
      /// @param message The message to send
      void sendToParent(const QString & message);
      
      private slots:
      
      /// @brief Slot called when the simulation has finished.
      void simulationFinished();
      
      /// @brief Slot called when a new frame comes from a MPI intercommunicator.
      
      /// @param senderComm Sender MPI communicator
      /// @param frameInfo Received frame information
      void newFrame(const MPI::Intercomm & senderComm, 
                    const COOLFluiD::Config::BuilderParserFrameInfo & frameInfo);
      
    private:
      
      /// @brief Manager-Worker m_protocol rules
      COOLFluiD::Config::ManagerWorkerProtocol m_protocol;
      
      /// @brief MPI listener
      MPIListener m_listener;
      
      /// @brief Server simulation
      ServerSimulation * m_srvSimulation;
      
      /// @brief Subsystem name
      QString m_subSystemName;
      
      /// @brief Subsystem type
      QString m_subSystemType;
      
      /// @brief Manager communicator
      MPI::Intercomm m_managerComm;
      
      /// @brief Communicators to other subsystems
      
      /// The key is the subsystem name. The value is the associated communicator.
      QHash<QString, MPI::Intercomm> m_comms;
      
      /// @brief The process rank.
      int m_rank;
      
      /// @brief Builds and sends a frame to the manager.
      
      /// @param fi Frame information
      void buildAndSend(const COOLFluiD::Config::BuilderParserFrameInfo & fi);
      
      /// @brief Builds and send an ACK frame
      
      /// @param type Frame type to ACK.
      void ack(COOLFluiD::Config::ManagerWorkerFrameType type);
      
      /// @brief Runs the simulation
      void runSimulation();
      
      /// @brief Opens a port
      
      /// @param remoteRole Remote subsystem name
      void openPort(const std::string & remoteRole);
      
      /// @brief Connects to a port
      
      /// @param portName Port name to connect to.
      /// @param remoteRole Remote subsystem name
      void connectToPort(const std::string & portName, const std::string & remoteRole);
      
      void setStatus(COOLFluiD::Common::WorkerStatus::Type status);
      
      void setBarrier(MPI::Intercomm comm);
      
      void setBarrier(MPI::Intercomm comm, COOLFluiD::Common::WorkerStatus::Type newStatus);
      
      void createSimulator();
      
      
    }; // class SimulationWorker
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_SimulationWorker_h
