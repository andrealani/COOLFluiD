#ifndef COOLFluid_server_ServerNetworkComm_h
#define COOLFluid_server_ServerNetworkComm_h

/////////////////////////////////////////////////////////////////////////////

#include <QObject>
#include <QAbstractSocket>
#include <QDomDocument>
#include <QList>
#include <QMutex>

#include "Config/BuilderParserFrameInfo.hh"
#include "ClientServer/network/NetworkProtocol.hh"
#include "ClientServer/network/HostInfos.hh"

class QHostAdress;
class QTcpServer;
class QTcpSocket;
class QString;

namespace COOLFluiD
{
  namespace network
  {
    class NetworkException;
    struct FrameInfos;
  }
  
  namespace server
  {
    
    class ServerSimulation;
    
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief This class is the server network level.
    
    /// For all signals emitted by this class, the first parameter identifies the
    /// client that caused the signal to be emitted. A same client will always
    /// have the same id. When client disconnects, its id will never be given
    /// to another client.
    
    /// @author Quentin Gasper. 
    class ServerNetworkComm : public QObject
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor.
      
      /// Creates a new m_socket with the address and port provided. The m_socket
      /// @c client is set to @c NULL.
      /// @param hostAddress Server address.
      /// @param port Socket port.
      /// @throw NetworkException Throws a NetworkException if the server cannot 
      /// listen to the given address/port.
      ServerNetworkComm(QString hostAddress = "127.0.0.1", quint16 port = 62784);
      
      /// @brief Destructor.
      
      /// Closes the sockets before the object is deleted.
      ~ServerNetworkComm();
      
      /// @brief Gives the number of bytes recieved.
      
      /// @return Returns the number of bytes recieved.
      int getBytesRecieved() const;
      
      /// @brief Gives the number of bytes sent.
      
      /// @return Returns the number of bytes sent.
      int getBytesSent() const;
      
      /// @brief Sends directory contents to a client
      
      /// @param clientId Client ID. If -1, the frame will be send to all client.
      /// @param path Directory path
      /// @param dirs List of sub-directories. May be empty.
      /// @param files List of files. May be empty.
      /// @return Returns @c true if the frame was built and sent successfully; 
      /// otherwise returns @c false.
      bool sendDirContents(int clientId, const QString & path,
                           const QStringList & dirs, const QStringList & files);
      
      /// @brief Send abstract types to a client.
      
      /// @param clientId Client ID. If -1, the frame will be send to all client.
      /// @param typeName Type name.
      /// @param typeList List of abstract types.
      /// @return Returns @c true if the frame was built and sent successfully; 
      /// otherwise returns @c false.
      bool sendAbstractTypes(int clientId, const QString & typeName,
                             const QStringList & typeList);
      
      /// @brief Send concrete types to a client.
      
      /// @param clientId Client ID. If -1, the frame will be send to all client.
      /// @param typeName Type name.
      /// @param typeList List of concrete types.
      /// @return Returns @c true if the frame was built and sent successfully; 
      /// otherwise returns @c false.
      bool sendConcreteTypes(int clientId, const QString & typeName,
                             const QStringList & typeList);
      
      /// @brief Send the host list to a client
      
      /// @param clientId Client ID. If -1, the frame will be send to all client.
      /// @param hosts Host list to send
      /// @return Returns @c true if the frame was built and sent successfully; 
      /// otherwise returns @c false.
      bool sendHostList(int clientId, 
                        const QList<COOLFluiD::network::HostInfos> & hosts);
      
      /// @brief Sends the tree to a client.
      
      /// @param clientId Client id, or -1 to send to all clients.
      /// @param tree The tree.
      /// @throws InvalidValueException If @c clientId does not exist.
      /// @throw UnknownClientIdException if Client id is unknown.
      void sendTree(int clientId, const QDomDocument & tree);
      
      /// @brief Sends the tree to a client.
      
      /// @param clientId Client id, or -1 to send to all clients.
      /// @param tree The tree.
      /// @throws InvalidValueException If @c clientId does not exist.
      /// @throw UnknownClientIdException if Client id is unknown.
      void sendTree(int clientId, const XMLNode & tree);
      
      /// @brief Sends an ACK/NACK.
      
      /// @param clientId Client id, or -1 to send to all clients.
      /// @param success If @c true an "ACK" frame is built. Otherwise, a "NACK" 
      /// frame is built.
      /// @param type Type of the frame to ACK/NACK.
      /// @throw UnknownClientIdException if Client id is unknown.
      void sendAck(int clientId, bool success, COOLFluiD::network::NetworkFrameType type);
      
      /// @brief Sends an already built frame.
      
      /// @param clientId Client id, or -1 to send to all clients.
      /// @param frame Frame to send
      /// @throw UnknownClientIdException if Client id is unknown.
      void send(int clientId, const QString & frame) ;
      
      /// @brief Sends an error message to a client 
      
      /// @param clientId Client id, or -1 to send to all clients.
      /// @param message Message to send
      /// @throw UnknownClientIdException if Client id is unknown.
      void sendError(int clientId, const QString & message);
      
      /// @brief Sends a message to a client
      
      /// @param clientId Client id, or -1 to send to all clients.
      /// @param message Message to send
      /// @throw UnknownClientIdException if Client id is unknown.
      void sendMessage(int clientId, const QString & message);
      
      void sendStatus(int clientId, const QString & subSysName, int rank, 
                      const QString & status);
      
      void sendSubSystemList(int clientId, const QStringList & subSystemList);
      
      private slots :
      
      /// @brief Slot called when a new client tries to connect.
      void newClient();
      
      /// @brief Slot called when new data are available on one of the 
      /// client sockets.
      void newData();
      
      /// @brief Slot called when the client has been disconnected.
      void clientDisconnected();
      
      /// @brief Slot called when a message (i.e. output forwarding) comes 
      /// from the simulator.
      
      /// Sends this message to all clients.
      /// @param message The message.
      void message(const QString & message);
      
      /// @brief Slot called when an error comes from the simulator.
      
      /// Sends this error message to all clients.
      /// @param message The error message.
      void error(const QString & message);
      
    signals:
      
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
      void renameNode(int clientId, const QString & nodePath, const QString & newName);
      
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
      /// @param clientId Client that requested the simulation to be run. All 
      /// errors message will be sent to this client. If the simulation has been
      /// successfully started, an ack for network frame type 
      /// @c TYPE_RUN_SIMULATION is sent to all clients. This parameter may 
      /// be @c NULL.
      void runSimulation(int clientId);
      
      /// @brief Signal emitted when the client want to activate the simulation
      
      /// @param clientId Client ID.
      /// @param nbProcs Number of workers to spawn
      /// @param hosts Hosts on which the workers will be spawned
      void activateSimulation(int clientId, unsigned int nbProcs, const QString & hosts);
      
      /// @brief Signal emitted when the client want to deactivate the simulation
      
      /// @param clientId Client ID.
      void deactivateSimulation(int clientId);
      
      /// @brief Shuts the server down
      
      /// Exits the Qt's event loop.
      /// @param clientId Client id, or -1 to send notifactions to all clients.
      void shutdownServer(int clientId);
      
      /// @brief Signal emitted when the client request the host list
      
      /// @param clientId Client ID.
      void getHostList(int clientId);
      
      void getSubSysList(int clientId);
      
    private:
      
      /// @brief The server m_socket.
      
      /// Used to accept connections.
      QTcpServer * m_server;
      
      /// @brief The server m_socket for the local loop.
      
      /// Used to accept connections coming from "localhost" (local loop).
      QTcpServer * m_localSocket;
      
      /// @brief Size of the frame that is being read.
      
      /// If the value is 0, no frame is currently being recieved.
      quint32 m_blockSize;
      
      /// @brief The client sockets.
      
      /// The key is pointer to the m_socket. ...
      QHash<QTcpSocket *, QDomNode> m_clients;
      
      /// @brief Number of bytes recieved.
      int m_bytesRecieved;
      
      /// @brief Number of bytes sent.
      int m_bytesSent;
      
      /// @brief Last ID given to a client
      int m_lastClientId;
      
      /// @brief Hash map that associated an ID to a client m_socket
      
      /// The key is the client ID. The value is the client m_socket.
      QHash<int, QTcpSocket *> m_clientIds;
      
      /// @brief Network m_protocol rules
      COOLFluiD::network::NetworkProtocol m_networkProtocol;
      
      /// @brief Builds and sends a frame to a client
      
      /// @param client Client to which the frame has to be sent. If @c NULL, the
      /// frame is sent to all clients.
      /// @param frameInfos Frame information
      /// @return Returns @c true if the frame was built and sent successfully; 
      /// otherwise returns @c false.
      bool buildAndSend(QTcpSocket * client, 
                        const COOLFluiD::Config::BuilderParserFrameInfo & frameInfos);
      
      /// @brief Sends a message to a client.
      
      /// @param client Client m_socket to use. If @c NULL, the message will be 
      /// sent to all clients.
      /// @param message Message to send.
      /// @return Returns @c true if the frame was built and sent successfully; 
      /// otherwise returns @c false.
      bool sendMessage(QTcpSocket * client, const QString & message);
      
      /// @brief Sends an error message to a client.
      
      /// @param client Client m_socket to use. If @c NULL, the error message will 
      /// be sent to all clients.
      /// @param message Error message to send.
      void sendError(QTcpSocket * client, const QString & message);
      
      /// @brief Sends a message to a client.
      
      /// @param client Client m_socket to use. If @c NULL, the frame will be sent 
      /// to all clients.
      /// @param frame Frame to send.
      /// @return Returns the number of bytes sent.
      int send(QTcpSocket * client, const QString & frame);
      
      /// @brief Sends an ACK/NACK.
      
      /// @param client Client m_socket to use. If @c NULL, the frame will be sent 
      /// to all clients.
      /// @param success If @c true an "ACK" frame is built. Otherwise, a "NACK" 
      /// frame is built.
      /// @param type Type of the frame to ACK/NACK.
      void sendAck(QTcpSocket * client, bool success, COOLFluiD::network::NetworkFrameType type);
      
      /// @brief Retrieves a client m_socket from its id.
      
      /// @param clientId Client id
      /// @return Returns a pointer to the m_socket, or @c NULL if client id was
      /// -1 (all clients).
      /// @throw UnknownClientIdException if Client id is unknown.
      QTcpSocket * getSocket(int clientId) const;
      
    }; // class ServerNetworkComm
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace network
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluid_server_ServerNetworkComm_h
