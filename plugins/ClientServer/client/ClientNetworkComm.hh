#ifndef COOLFluid_ClientServer_ClientNetworkComm_h
#define COOLFluid_ClientServer_ClientNetworkComm_h

/////////////////////////////////////////////////////////////////////////////

#include <QObject>
#include <QAbstractSocket>
#include <QStringList>

#include "Config/BuilderParserFrameInfo.hh"
#include "ClientServer/network/HostInfos.hh"
#include "ClientServer/network/NetworkException.hh"
#include "ClientServer/network/NetworkProtocol.hh"

class QDomDocument;
class QDomNode;
class QString;
class QTcpServer;
class QTcpSocket;

namespace COOLFluiD
{
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class represents the client network level.
    
    /// It operates mainly using Qt slots/signals system. Each time a frame 
    /// arrives through the m_socket, the appropriate signal is thrown. Frames 
    /// to send are built using the BuilderParser system.
    
    /// @author Quentin Gasper.
    
    class ClientNetworkComm : public QObject
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor.
      
      /// The m_socket @c #commClient is set to @c NULL.
      ClientNetworkComm();
      
      /// @brief Destructor.
      
      /// Closes the sockets and free all allocated memory before the object 
      /// is deleted.
      ~ClientNetworkComm();
      
      /// @brief Attempts to connect the client to the server.
      
      /// When this method returns, the m_socket is not open yet. The signal
      /// @c connected() will be emitted when the first frame 
      /// arrives.
      /// @param hostAddress Server address.
      /// @param port Socket port number.
      /// @param m_skipRefused Value of @c #m_skipRefused during the attempt.
      void connectToServer(const QString & hostAddress = "127.0.0.1",
                           quint16 port = 62784, bool skipRefused = false);
      
      /// @brief Disconnects from the server, then closes.
      
      /// After calling this method, @c #m_requestDisc is @c true.
      /// @param shutServer If @c true, a request to shut down the server is sent.
      void disconnectFromServer(bool shutServer);
      
      /// @brief Sends a request to the server to get the tree.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendGetTree();
      
      /// @brief Sends a request to the server to get the host list.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendGetHostList();
      
      /// @brief Sends a request to the server to add a node.
      
      /// The node parents indicate the path in the tree and all parents must 
      /// already exist in the tree on the server, otherwise the server will 
      /// send back an error.
      /// @param node Node to add.
      /// @param type Concrete type of the node.
      /// @param absType Abstract type of the node.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendAddNode(const QDomNode & node, const QString & type, 
                       const QString & absType);
      
      /// @brief Sends a request to the server to modify a node.
      
      /// @c data can be built using @c #TreeModel::buildModification() method.
      /// @param data Data about the modification.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendModifyNode(const QDomDocument & data);
      
      /// @brief Sends a request to the server to delete a node.
      
      /// The node parents indicate the path in the tree and all parents and 
      /// node to delete must exist in the tree on the server, otherwise the 
      /// server will send back an error.
      /// @param node Node to delete.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendDeleteNode(const QDomNode & node);
      
      /// @brief Sends a request to the server to close the file.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendCloseFile();
      
      /// @brief Sends a request to the server to rename a node.
      
      /// The node parents indicate the path in the tree and all parents and 
      /// node to rename must exist in the tree on the server, otherwise the 
      /// server will send back an error. The server will also send back an 
      /// error if the another node with the same name as the new already 
      /// exists. If the node name and the new name are the same, there is no 
      /// error.
      /// @param node Node to rename.
      /// @param newName Node new name.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendRenameNode(const QDomNode & node, const QString & newName);
      
      /// @brief Sends a request to the server to get the abstract types for a 
      /// specified type.
      
      /// @param typeName Type name
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendGetAbstractTypes(const QString & typeName);
      
      /// @brief Sends a request to the server to get the concrete types for a 
      /// specified abstract type.
      
      /// @param typeName Name of the abstract type one want to get the 
      /// concrete types list.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendGetConcreteTypes(const QString & typeName);
      
      /// @brief Sends a request to the server to open a case file.
      
      /// @param filename File to open
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendOpenFile(const QString & filename);
      
      /// @brief Sends a request to the server to run a simulation.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendRunSimulation();
      
      /// @brief Sends a request to open a directory and read its content.
      
      /// @param dirname Directory name to open.
      /// @param includeFiles If @c true, sub-directories and files will be 
      /// requested. If @c false, only sub-directories will be requested.
      /// @param extensions List of wanted extensions. According to the network
      /// m_protocol, this list may be empty.
      /// @param includeNoExtension If @c true, files without any extension will 
      /// be requested. If @c false, they will not.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendOpenDir(const QString & dirname,
                       bool includeFiles,
                       const QStringList & extensions = QStringList(),
                       bool includeNoExtension = false);
      
      /// @brief Indicates wether a connection is established to a server.
      
      /// @return Returns @c true if the client is connected to a server; 
      /// otherwise, returns @c false.
      bool isConnected() const;
      
      /// @brief Sends a configuration to the server.
      
      /// @param document Configuration to send.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendConfig(const QDomDocument & document);
      
      /// @brief Sends a request to the server to create a new directory.
      
      /// @param path Path to the parent directory.
      /// @param name New directory name.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendCreateDir(const QString & path, const QString & name);
      
      /// @brief Sends a request to the server to save a configuration tree.
      
      /// @param path File path to write.
      /// @param config Configuration tree to save
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendSaveConfig(const QString & path, const QDomDocument & config);
      
      /// @brief Sends a request to the server to activate the simulation.
      
      /// @param nbProcs Number of workers to spawn.
      /// @param hosts Host list on which the workers will be spawn. May be
      /// empty.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendActivateSimulation(int nbProcs, const QString & hosts);
      
      /// @brief Sends a request to the server to deactivate the simulation.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool sendDeactivateSimulation();
      
      bool sendGetSubSystemList();
      
      public slots :
      
      /// @brief Slot called when there is an error on the m_socket.
      void newData();
      
      /// @brief Slot called when the connection has been successfully established.
      void connectionEstablished();
      
      /// @brief Slot called when the connection has been broken.
      void disconnected();
      
      /// @brief Slot called when there is an error on the m_socket.
      
      /// @param err Error that occured.
      void socketError(QAbstractSocket::SocketError err);
      
    signals:
      
      /// @brief Signal emitted when the server sends the tree.
      
      /// @param document The tree.
      void newTree(const QDomDocument & document);
      
      /// @brief Signal emitted when the m_socket has been closed due to a 
      /// network error.
      
      /// The signal is not emitted if the user resquested a disconnection (if 
      /// @c #m_requestDisc is @c true ).
      void disconnectFromServer();
      
      /// @brief Signal emitted when a connection has been successfully 
      /// established between the client and the server. 
      
      /// The signal is emitted exactly once when the first frame is 
      /// recieved from the server.
      void connected();
      
      /// @brief Signal emitted when the server sends an abstract types list 
      /// of a concrete type.
      
      /// @param types Abstract types list. Each element is a type.
      void abstractTypes(const QStringList & types);
      
      /// @brief Signal emitted when the server sends an concrete types list.
      
      /// @param types Concrete types list. Each element is a type.
      void concreteTypes(const QStringList & types);
      
      /// @brief Signal emitted when the server sends an ACK (acknowledgement) 
      /// for a specified type of frame.
      
      /// @param type Type of the acknowledged frame.
      void ack(COOLFluiD::network::NetworkFrameType type);
      
      /// @brief Signal emitted when the server sends an NACK 
      /// (non-acknowledgement) for a specified type of frame.
      
      /// @param type Type of the non-acknowledged frame.
      void nack(COOLFluiD::network::NetworkFrameType type);
      
      /// @brief Signal emitted when the server sends a directory contents.
      
      /// @param path Absolute path of the directoy of which contents belong to.
      /// @param dirs Directories list. Each element is a directory.
      /// @param files Files list. Each element is a file.
      void dirContents(const QString & path, const QStringList & dirs, 
                       const QStringList & files);
      
      /// @brief Signal emitted when the server sends the host list.
      
      /// @param infos List of hosts.
      void hostList(const QList<COOLFluiD::network::HostInfos> & infos);
      
      void simulationStatus(const QString & subSysName, int rank, 
                            const QString & status);
      
      void subsystemList(const QStringList & subSystems);
      
    private: // methods
      
      /// @brief Sends a frame to the server.
      
      /// All @e sendXXX methods of this class call this method to 
      /// send their frames.
      /// @param frame Frame to send.
      int send(const QString & frame) const;
      
      /// @brief Builds and sends a frame to the server.
      
      /// The frame is sent only if build was successful; otherwise @c #error()
      /// signal is emitted with error message.
      /// @param frameInfos Frame information
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool buildAndSend(const COOLFluiD::Config::BuilderParserFrameInfo & frameInfos);
      
      /// @brief Builds and sends a frame to the server.
      
      /// This is an overloaded method, provided for convinience. It can be used 
      /// to send a frame that does not need neither attributes nor data.
      /// @param type Frame type.
      /// @return Return @c true if the frame was successfully built and sent; 
      /// otherwise returns @c false.
      bool buildAndSend(COOLFluiD::network::NetworkFrameType type);
      
      /// @brief Recursive method that builds a @c QSting with 
      /// all parent nodes names of a given node. 
      
      /// The path is UNIX path-like formatted.
      /// @param node Node from which the parents are returned
      /// @return Returns the built string
      QString getNodePath(const QDomNode & node) const;
      
      bool checkConnected();
      
    private: // data
      
      /// @brief Socket used to communicate with the server.
      QTcpSocket * m_socket;
      
      /// @brief Size of the frame that is being read.
      
      /// If the value is 0, no frame is currently being recieved.
      quint32 m_blockSize;
      
      /// @brief Indicates wether the upper level requested a disconnection.
      bool m_requestDisc;
      
      /// @brief Indicates wether the m_socket is open and connected to the 
      /// server.
      bool m_connectedToServer;
      
      /// @brief Indicates wether a "Connection refused" error must be skip.
      
      /// If @c true when a "Connection refused" error occurs, it is 
      /// skipped and this attribute is set to @c false.
      bool m_skipRefused;
      
      /// @brief Network m_protocol rules
      COOLFluiD::network::NetworkProtocol m_protocol;
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namspace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluid_ClientServer_ClientNetworkComm_h
