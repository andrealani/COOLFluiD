#ifndef COOLFluiD_server_ServerKernel_h
#define COOLFluiD_server_ServerKernel_h

/////////////////////////////////////////////////////////////////////////////

#include <QObject>
#include <QThread>
#include <QList>
#include <QDomDocument>
#include <mpi.h>

#include "Common/PE.hh"
#include "ClientServer/network/HostInfos.hh"
#include "ClientServer/server/ServerSimulation.hh"
#include "ClientServer/server/SimulationManager.hh"

namespace COOLFluiD
{
  namespace server
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    class ServerNetworkComm;
    class ServerSimulation;
    class SimulationManager;
    
    /// @brief Manages everything that is not related to the network layer
    /// or the simulation management.
    
    /// @author Quentin Gasper
    
    class ServerKernel : public QObject
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor
      
      /// Builds the network layer
      /// @param hostame Host name
      /// @param portNumber Port number
      ServerKernel(const QString & hostname, quint16 portNumber, 
                   const QList<COOLFluiD::network::HostInfos> & hostList);
      
      /// @brief Destructor
      
      /// Free all allocated memory.
      ~ServerKernel();
      
      private slots:
      
      /// @brief Slot called when a new client connects
      
      /// Sends server status (file open, simulation running) to the new client.
      /// @param clientId New client id.
      void newClient(int clientId);
      
      /// @brief Sends the tree to a client.
      
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      /// @param tree The tree.
      void getTree(int clientId);
      
      /// @brief Configures the simulator
      
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      /// @param config The configuration to give to the simulator
      void configureSimulator(int clientId, const QDomDocument & config);
      
      /// @brief Requests the simulator to modify a node
      
      /// The simulator must have been configured.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      /// @param data Data about the modification
      void modifyNode(int clientId, const QDomDocument & data);
      
      /// @brief Requests the simulator to delete a node
      
      /// The simulator must have been configured and a simulation can not be 
      /// running.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      /// @param data Data about the modification
      void deleteNode(int clientId, const QString & nodePath);
      
      /// @brief Requests the simulator to add a node
      
      /// The simulator must have been configured and a simulation can not be 
      /// running.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      /// @param parentPath Path of the new node parent 
      /// @param name Name of the new node
      /// @param type Node type
      /// @param absType Node abstract type
      void addNode(int clientId, const QString & parentPath, 
                   const QString & name, const QString & type, 
                   const QString & absType);
      
      /// @brief Requests the simulator to rename a node
      
      /// The simulator must have been configured and a simulation can not be 
      /// running.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      /// @param m_nodePath Node path
      /// @param newName Node new name
      void renameNode(int clientId, const QString & nodePath, 
                      const QString & newName);
      
      /// @brief Sends the abstract types list.
      
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      
      /// @param typeName Type name of which the abstract types are requested.
      void getAbstractTypes(int clientId, const QString & typeName);
      
      /// @brief Sends the concrete types list.
      
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      
      /// @param typeName Abstract type name of which the abstract types are
      /// requested.
      /// @param typesList Concrete types list.
      void getConcreteTypes(int clientId, const QString & typeName);
      
      /// @brief Creates a new directory
      
      /// Creates all the parent that does not exist.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      
      /// @param dirPath Directory path
      /// @param name Directory name
      void createDirectory(int clientId, const QString & dirPath, 
                           const QString & name);
      
      /// @brief Reads a directory contents.
      
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      
      /// @param dirPath Path to the directory to read.
      /// @param extensions List of wanted extensions. According to the network
      /// m_protocol, this list may be empty.
      /// @param includeFiles If @c true, sub-directories and files returned. 
      /// If @c false, only sub-directories will be returned.
      /// @param includeNoExtension If @c true, files without any extension will 
      /// be returned. If @c false, they will not.
      void openDirectory(int clientId, const QString & dirPath, 
                         const QStringList & extensions, bool includeFiles,
                         bool includeNoExtension);
      
      /// @brief Saves a configuration to a file
      
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      
      /// @param filename File to write
      /// @param config Configuration
      /// @note If the file already exists, it is overwritten.
      void saveConfiguration(int clientId, const QString & filename, 
                             const QDomDocument & config);
      
      /// @brief Requests to the simulator to open a file.
      
      /// This method returns when the file is successfully opened or when an 
      /// error has occured.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      
      /// @param filename Name of the file to open.
      void openFile(int clientId, const QString & filename);
      
      /// @brief Closes the file
      
      /// This destroys the simulator and create a new one.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      void closeFile(int clientId);
      
      /// @brief Requests to the simulator to run the simulation.
      
      /// Starts the simulator thread and returns immediately.
      /// @param client Client that requested the simulation to be run. All 
      /// errors message will be sent to this client. If the simulation has been
      /// successfully started, an ack for network frame type 
      /// @c TYPE_RUN_SIMULATION is sent to all clients. This parameter may 
      /// be @c NULL.
      void runSimulation(int clientId);
      
      /// @brief Shuts the server down
      
      /// Exits the Qt's event loop.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      void shutdownServer(int clientId);
      
      /// @brief Forwards a message from the simulator to the network layer
      /// @param message Message to send
      void message(const QString & message);
      
      /// @brief Forwards an error message from the simulator to the network layer
      /// @param message Error message to send
      void error(const QString & message);
      
      /// @brief Slot called when the simulation is done.
      void simulationFinished();
      
      /// @brief Slot called when the client wants the host list
      
      /// @param clientId Client ID.
      void getHostList(int clientId);
      
      /// @brief Slot called when the client wants to activate the simulation
      
      /// @param clientId Client ID.
      /// @param nbProcs Number of workers to spawn
      /// @param hosts Hosts on which the workers will be spawned.
      void activateSimulation(int clientId, unsigned int nbProcs, 
                              const QString & hosts);
      
      /// @brief Slot called when the client wants to deactivate the simulation
      /// @param clientId Client ID.
      void deactivateSimulation(int clientId);
      
      void simulationStatus(const QString & subSysName, int rank, 
                            const QString & status);
      
      void getSubSysList(int clientId);
      
      void simulationTree(const XMLNode & tree);
      
      void spawned();
      
    private:
      /// @brief The default path for the file browsing.
      
      /// The default path is the current directory (./).
      const QString DEFAULT_PATH;
      
      /// @brief The simulation.
      ServerSimulation * m_srvSimulation;
      
      /// @brief The network communication
      ServerNetworkComm * m_commServer;
      
      /// @brief XML tree with available types, their concrete and abstract 
      /// sub-types.
      
      /// The content of this tree is read from "./TypesList.xml" file
      QDomDocument m_types;
      
      /// @brief Indicates wether a file is already open.
      
      /// If @c true, a file is already open.
      bool m_fileOpen;
      
      /// @brief Indicates wether the simulation is running.
      
      /// If @c true, the simulation is running.
      bool m_simRunning;
      
      bool m_active;
      
      /// @brief Simulation manager
      SimulationManager m_simulationManager;
      
      /// @brief Simulation tree
      /// @todo this attribute should be removed
      XMLNode m_simTree;
      
      /// @brief Host list
      QList<COOLFluiD::network::HostInfos> m_hostList;
      
      /// @brief Creates a new simulator with the name @c name.
      
      /// If a simulator already exists, it is not destroyed. If a file is open, 
      /// it is not closed.
      /// @param name Name of the new simulator.
      void createSimulator(const QString & name);
      
      /// @brief Reads a directory contents.
      
      /// @param directory Directory to read.
      /// @param extensions List of wanted extensions. According to the network
      /// m_protocol, this list may be empty.
      /// @param includeFiles If @c true, sub-directories and files returned. 
      /// If @c false, only sub-directories will be returned.
      /// @param includeNoExtension If @c true, files without any extension will 
      /// be returned. If @c false, they will not.
      /// @param dirsList Reference of a @c QStringList where sub-directories 
      /// names will be stored.
      /// @param filesList Reference of a @c QStringList where files names will 
      /// be stored.
      /// @return Returns @c true if the directory has been correctly read. 
      /// Otherwise, returns @c false (@c dirsList and @c filesList are not 
      /// modified in this case).
      bool getDirContent(const QString & directory, 
                         const QStringList & extensions,
                         bool includeFiles,
                         bool includeNoExtension,
                         QStringList & dirsList, 
                         QStringList & filesList) const;
      
      /// @brief Connect a signal from network layer to slots
      
      /// @param signal The signal
      /// @param slot The slot
      /// @return Returns @c true if the signal was well connected to the slot; 
      /// otherwise, returns @c false.
      bool connSig(const char * signal, const char * slot);
      
      void setStatus(COOLFluiD::Common::WorkerStatus::Type status);
      
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_ServerKernel_h
