#ifndef COOLFluiD_client_ClientKernel_h
#define COOLFluiD_client_ClientKernel_h

#include <QMap>
#include <QObject>

#include "ClientServer/network/NetworkFrameType.hh"
#include "ClientServer/network/HostInfos.hh"
#include "ClientServer/treeview/TreeModel.hh"

class QProcess;
class QString;
class QTimer;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    class StatusModel;
    
    //////////////////////////////////////////////////////////////////////////////
    
    class ClientNetworkComm;
    
    class ClientKernel : public QObject
    {
      Q_OBJECT
      
    public:
      
      static ClientKernel * getInstance();
      
      void setTreeModel(COOLFluiD::treeview::TreeModel * treeModel);
      
      void setStatusModel(StatusModel * statusModel);
      
      COOLFluiD::treeview::TreeModel * getTreeModel() const;
      
      StatusModel * getStatusModel() const;
      
      void readDir(const QModelIndex & index, const QString & dirname,
                   bool includeFiles, 
                   const QStringList & extensions = QStringList(),
                   bool includeNoExtension = false);
      
      void createDir(const QModelIndex & index, const QString & parentPath, 
                     const QString & dirName);
      
      void openFile(const QModelIndex & index, const QString & filename);
      
      private slots:
      
      /// @brief Sets a new TreeModel to the the treeview. 
      
      /// The model is built from the given document.
      /// @param domDocument Document to use to build the new TreeView.
      void newTree(const QDomDocument & m_domDocument);
      
      /// @brief Slot called when the client is connected to the server.
      void connected();
      
      /// @brief Slot called when the network level recieves abstract types 
      /// list.
      
      /// @param types Recieved abstract types list.
      void abstractTypes(const QStringList & types);
      
      /// @brief Slot called when the network level recieves concrete types 
      /// list.
      
      /// @param types Recieved concrete types list.
      void concreteTypes(const QStringList & types);
      
      /// @brief Tries to connect to the server.
      
      /// During the waiting for the server to launch through an SSH 
      /// connection, this slot is called at every timeout of
      /// @c #timer and tries to connect to the server.
      void tryToConnect();
      
      /// @brief Slot called when server sends host list
      
      /// @param infos The list
      void hostList(const QList<COOLFluiD::network::HostInfos> & infos);
      
      void subsystemList(const QStringList & subSystems);
      
      /// @brief Slot called when the user wants to disconnect from the server.
      
      /// The user is invited to choose between shutdown the server, just 
      /// disconnect from it or cancel the action. If the user don't select 
      /// "Cancel", this slot destroys @c #communication. Both pointers are set 
      /// to @c NULL. This slot destroys the @c ClientNetworkComm object and the model.
      void disconnectFromServer();
      
      /// @brief Signal emitted when an ACK (acknowledgement) arrives for a frame.
      
      /// @param type Type of the acknowledged frame (conforming to type 
      /// defined by FrameBuilder class).
      void ack(COOLFluiD::network::NetworkFrameType type);
      
      /// @brief Signal emitted when a NACK (non-acknowledgement) arrives for a
      /// frame.
      
      /// @param type Type of the non-acknowledged frame (conforming to type
      /// defined by FrameBuilder class).
      void nack(COOLFluiD::network::NetworkFrameType type);
      
      void simulationStatus(const QString & subSysName, int rank, 
                            const QString & status);
      
      void sshError();
      
      public slots:
      
      /// @brief Slot called when user wants to connect a simulation to its server.
      
      /// @param index Simulation index
      /// @param info Connection information
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void connectSimulation(const QModelIndex & index, 
                             const COOLFluiD::treeview::TSshInformation & info);
      
      /// @brief Slot called when user wants to disconnect a simulation from 
      /// its server.
      
      /// @param index Simulation index
      /// @param shutServer If @c true, the server shuts down. Otherwise, it keeps
      /// running after the client has been disconnected.
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void disconnectSimulation(const QModelIndex & index, bool shutServer);
      
      /// @brief Slot called when user wants to run a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void runSimulation(const QModelIndex & index);
      
      /// @brief Slot called when user wants to stop a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void stopSimulation(const QModelIndex & index);
      
      /// @brief Slot called when user wants to activate a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void activateSimulation(const QModelIndex & index);
      
      /// @brief Slot called when user wants to deactivate a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void deactivateSimulation(const QModelIndex & index);
      
      /// @brief Slot called when the user wants to add a child node of a
      /// selected abstract type.
      
      /// This slot sends a request to the server to get the concrete types 
      /// corresponding to the selected abstract type. It returns without 
      /// waiting an m_answer.
      void addNode(const QString & abstractType);
      
      /// @brief Slot called when the user wants to rename an object.
      
      /// A request is sent to the server if and only if the new name is not 
      /// empty and it is different to the old one.
      
      /// @param node Node to rename
      /// @param newName Node new name
      void renameNode(const QDomNode & node, const QString & newName);
      
      /// @brief Slot called when the user wants to delete an object.
      
      /// The slot sends the request to the server without asking the user to
      /// confirm (should be fixed in the future).
      
      /// @param node Node to delete
      void deleteNode(const QDomNode & node);
      
      void commitChanges(const QDomDocument & doc);
      
    signals:
      
      /// @brief Signal emitted when the server sends a directory contents.
      
      /// @param path Absolute path of the directoy of which contents belong to.
      /// @param dirs Directories list. Each element is a directory.
      /// @param files Files list. Each element is a file.
      void dirContents(const QString & path, const QStringList & dirs, 
                       const QStringList & files);
      
      void acked(COOLFluiD::network::NetworkFrameType type);
      
    private: // methods
      
      ClientKernel();
      
      ~ClientKernel();
      
      /// @brief Retrieves a simulation index from a signal sender, if it is
      /// a @c QTcpSocket.
      
      /// @param sender The sender
      /// @return Returns the simulation index, or an invalid index if the
      /// simulation corresponding to the sender (or the sender was not a
      /// @c QTcpSocket ) could not be found.
      QModelIndex getSimIndex(QObject * sender) const;
      
      void connectToServer(const QModelIndex & simIndex);
      
      void launchServer(const QModelIndex & simIndex);
      
    private: // data
      
      static ClientKernel * m_instance;
      
      COOLFluiD::treeview::TreeModel * m_treeModel;
      
      StatusModel * m_statusModel;
      
      /// @brief Active network sockets.
      
      /// The key is the simulation index. The value is a pointer to the
      /// associated m_socket.
      /// @todo Replace QModelIndex by QPersistentModelIndex
      QMap<QModelIndex, ClientNetworkComm *> m_networkComms;
      
      QMap<QModelIndex, QTimer *> m_timers;
      
      QMap<QModelIndex, QProcess *> m_launchServerProcs;
      
      QMap<QModelIndex, COOLFluiD::treeview::TSshInformation> m_commsSshInfo;
      
    }; // class ClientKernel
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_ClientKernel_h
