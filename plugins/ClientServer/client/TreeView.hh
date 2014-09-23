#ifndef COOLFluiD_client_TreeView_h
#define COOLFluiD_client_TreeView_h


/////////////////////////////////////////////////////////////////////////////

#include <QTreeView>
#include <QList>
#include <QHash>

#include "ClientServer/treeview/TSshInformation.hh"

class QDomDocument;
class QDomNode;
class QMainWindow;
class QMenu;
class QMenuBar;
class QModelIndex;
class QSortFilterProxyModel;

namespace COOLFluiD
{
  namespace treeview { class TreeModel; }
  
  namespace client
  {
    class OptionPanel;
    
    /////////////////////////////////////////////////////////////////////////////
    
    enum TreeViewMenuActionType
    {
      ACTION_SIM_NEW_SIMULATION,
      ACTION_SIM_OPEN_SIMULATION,
      ACTION_SIM_END_SIMULATION,
      ACTION_SIM_CONNECT,
      ACTION_SIM_DISCONNECT,
      ACTION_SIM_UPDATE_TREE,
      ACTION_SIM_RUN_SIMULATION,
      ACTION_SIM_STOP_SIMULATION,
      ACTION_SIM_ACTIVATE_SIM,
      ACTION_SIM_DEACTIVATE_SIM,
      ACTION_OBJECT_ADD_NODE,
      ACTION_OBJECT_RENAME_NODE,
      ACTION_OBJECT_DELETE,
      ACTION_OBJECT_PROPERTIES
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class manages the tree m_view.
    
    /// @author Gasper Quentin
    
    class TreeView : public QTreeView
    {
      Q_OBJECT 
      
    public:
      /// @brief Constructor.
      
      /// @param optionsPanel Panel m_options of the selected node will be displayed.
      /// @param parent Parent window. May be @c NULL.
      /// @throws std::invalid_argument if @c optionsPanel is @c NULL.
      TreeView(OptionPanel * optionsPanel, QMainWindow * parent = NULL);
      
      /// @brief Destructor.
      
      /// Frees all allocated memory. Neither the m_options panel nor the parent 
      /// are destroyed.
      ~TreeView();
      
      /// @brief Sets abstract types list.
      
      /// @param abstractTypes Abstract types list.
      void setAbstractTypesList(const QStringList & abstractTypes);
      
      /// @brief Builds a new node.
      
      /// This method calls @link COOLFluiD::treeview::TreeModel::newChildToNode()
      /// TreeModel::newChildToNode() @endlink for the currently selected index.
      /// @param newNode New node name.
      /// @param doc The document the node will be added to. The presence of
      /// this parameter is due to the fact that a node can not exist if it 
      /// does not belong to a document.
      /// @return The built node.
      QDomNode newChildNode(const QString & newNode, QDomDocument & doc) const;
      
      /// @brief Changes the read-only mode.
      
      /// @param readOnly If @c true, the tree m_view will switch to the read-only 
      /// mode.
      void setReadOnly(bool readOnly);
      
      /// @brief Indicates whether the tree m_view is in read-only mode or not.
      
      /// @return Returns @c true if the tree m_view is in read-only mode. 
      /// Otherwise, returns @c false.
      bool isReadOnly() const;
      
      QAction * addSimToMenu(QMenu * menu);
      
      QAction * addSimToMenuBar(QMenuBar * menuBar);
      
      void setTreeModel(COOLFluiD::treeview::TreeModel * treeModel);
      
      COOLFluiD::treeview::TreeModel * getTreeModel() const;
      
    protected:
      
      /// @brief Method called when a mouse button is pressed in the treeview.
      
      /// This method overloads parent class method. Four cases are possible :
      /// @li If user right-clicks, a context menu is displayed. 
      /// @li If user left-clicks on another node than the currently selected one 
      /// @b and @c confirmChangeOptions() returns @c true, m_options in the 
      /// m_options panel are changed. 
      /// @li If user left-clicks on the selected node nothing is done. 
      /// @li Middle button has no effect.
      /// @param event Event that occured.
      virtual void mousePressEvent(QMouseEvent * event);
      
      virtual void keyPressEvent(QKeyEvent * event);
      
      private slots:
      
      /// @brief Slot called when user wants to add a node.
      void addNode();
      
      /// @brief Slot called when user wants to delete a node.
      void deleteNode();
      
      /// @brief Slot called when the user wants to rename an object.
      
      /// A request is sent to the server if and only if the new name is not 
      /// empty and it is different to the old one.
      void renameNode();
      
      /// @brief Slot called when the user wants to see an object properties.
      
      /// Properties are displayed in a message box.
      void showProperties();
      
      /// @brief Slot called when user wants to add an option to the currently 
      /// selected object.
      
      /// It is called when user clicks on any type in "Add an option" sub-menu of
      /// the context menu.
      void addOption();
      
      /// @brief Slot called user wants to commit modified m_options
      
      /// @param modOptions Modified m_options
      /// @param m_newOptions Newly added m_options
      void changesMade(const QDomDocument & modOptions, 
                       const QDomDocument & newOptions);
      
      /// @brief Slot called when user wants to create a new simulation.
      void newSimulation();
      
      /// @brief Slot called when user wants to load a case file.
      void openSimulation();
      
      /// @brief Slot called when user wants to delete a simulation.
      void endSimulation();
      
      /// @brief Slot called when user wants to connect a simulation to its server.
      void connectSimulation();
      
      /// @brief Slot called when user wants to disconnect a simulation from
      /// its server.
      void disconnectSimulation();
      
      /// @brief Slot called when user wants to update a simulation tree.
      void updateTree();
      
      /// @brief Slot called when user wants to run a simulation.
      void runSimulation();
      
      /// @brief Slot called when user wants to stop a running simulation.
      void stopSimulation();
      
      /// @brief Slot called when user wants to activate a simulation.
      void activateSimulation();
      
      /// @brief Slot called when user wants to deactivate simulation.
      void deactivateSimulation();
      
      void currentIndexChanged(const QModelIndex & index);
      
      void nodeActivated(const QModelIndex & index);
      
    signals:
      
      /// @brief Signal emitted when user wants to add a node.
      
      /// @param abstractType Abstract type of the new node.
      void addNode(const QString & abstractType);
      
      /// @brief Signal emitted when user wants to delete a node.
      
      /// @param node Node index in the model.
      void deleteNode(const QDomNode & node);
      
      /// @brief Signal emitted when user wants to rename a node.
      
      /// @param node Node index in the model.
      /// @param newName Node new name.
      void renameNode(const QDomNode & node, const QString & newName);
      
      /// @brief Signal emitted when user wants to commit modified m_options
      
      /// @param document Modification information
      void commitChanges(const QDomDocument & document);
      
      /// @brief Signal emitted when user wants to connect a simulation to its 
      /// server.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void connectSimulation(const QModelIndex & index, 
                             const COOLFluiD::treeview::TSshInformation & info);
      
      /// @brief Signal emitted when user wants to connect a simulation to its 
      /// server.
      
      /// @param index Simulation index
      /// @param shutServer If @c true, the server shuts down. Otherwise, it keeps
      /// running after the client has been disconnected.
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void disconnectSimulation(const QModelIndex & index, bool shutServer);
      
      /// @brief Signal emitted when user wants to load a case file.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void openSimulation(const QModelIndex & index);
      
      /// @brief Signal emitted when user wants to run a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void runSimulation(const QModelIndex & index);    
      
      /// @brief Signal emitted when user wants to stop a running simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void stopSimulation(const QModelIndex & index);
      
      /// @brief Signal emitted when user wants to activate a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void activateSimulation(const QModelIndex & index);
      
      /// @brief Signal emitted when user wants to deactivcate a simulation.
      
      /// @param index Simulation index
      /// @todo Replace QModelIndex by QPersistentModelIndex
      void deactivateSimulation(const QModelIndex & index);
      
    private:
      
      /// @brief Hashmap containing all available actions for menu m_items. 
      
      /// The key is a number defined by one of the constant integer attributes 
      /// of this class. The value is the action corresponding to this number.
      QHash<int, QAction *> m_actions;
      
      /// @brief List containing all actions for abstract types displayed in 
      /// the context menu. 
      
      /// These actions are not stored in @c #actions because they are not 
      /// identified by an integer and the list may be cleared several times 
      /// during application runtime.
      QList<QAction *> m_abstractTypesActions;
      
      /// @brief Simulation menu
      QMenu * m_simulationMenu;
      
      /// @brief Object menu
      QMenu * m_objectMenu;
      
      /// @brief Abstract types menu. 
      
      /// This is a sub-menu of the context menu.
      QMenu * m_mnuAbstractTypes;
      
      /// @brief "Add a child node" sub-menu.
      QMenu * m_mnuNewOption;
      
      /// @brief Panel used to display and modify m_options for a selected 
      /// object.
      OptionPanel * m_optionsPanel;
      
      /// @brief List of abstract types
      QStringList m_abstractTypes;
      
      /// @brief The model to display.
      COOLFluiD::treeview::TreeModel * m_treeModel;
      
      /// @brief Filter for the treeview. 
      
      /// Allows to switch between basic/advanced mode. The filter is used as the 
      /// treeview model. Its source is the tree model.
      QSortFilterProxyModel * m_modelFilter;
      
      /// @brief Indicates whether the tree m_view is in read-only mode or not.
      
      /// If @c true, the tree m_view is read-only mode. When it is read-only mode,
      /// all m_options in the context that may modify an object are disbaled. 
      /// "Delete", "Rename", "Add a child node" and "Add an option" m_items are then
      /// disabled.
      bool m_readOnly;
      
      /// @brief Buillds object menu
      void buildObjectMenu();
      
      /// @brief Builds the simulation menu
      void buildSimulationMenu();
      
      /// @brief Asks user to commit or rollback before changing m_options in 
      /// m_options panel.
      
      /// If modifications were committed, nothing is asked and the method
      /// immediately returns @c true. If the commit is requested by the user,
      /// it is processed by this method.
      /// @param index Node index on which user clicked. If it is equals to
      /// @c #currentIndexInPanel nothing is asked and the method
      /// immediately returns @c true.
      /// @param okIfSameIndex If @c false, the method checks if indexes are the 
      /// same. If @c true, no check is done. 
      /// @return Returns @c false if the user clicked on "Cancel" ; otherwise
      /// returns @c true.
      bool confirmChangeOptions(const QModelIndex & index, bool okIfSameIndex = false);
      
      /// @brief Saves the expanded or collapsed state of each node in a hash map.
      
      /// This is a recursive method. It calls itself on each child of the node 
      /// pointed by @c index. The key of the map is a string representing the 
      /// path of the node pointed by a model index as returned by 
      /// @c TreeModel::getNodePath(). The value of the map is a bool: @c true
      /// if the node is expanded or @c false if it is collapsed.
      /// @param index Node the save. All its children are saved too.
      /// @param map The map where the states will be saved to.
      /// @warning When calling the method for the first in the recursivity loop
      /// make sure that the map is empty to avoid map data inconsistency.
      void getChildrenNodesState(const QModelIndex & index, 
                                 QHash<QString, bool> & map);
      
      /// @brief Sets the expanded or collapsed state of each node in a hash map.
      
      /// This is a recursive method. It calls itself on each child of the node 
      /// pointed by @c index. The key of the map is a string representing the 
      /// path of the node pointed by a model index as returned by 
      /// @c TreeModel::getNodePath(). The value of the map is a bool: @c true
      /// if the node is expanded or @c false if it is collapsed.
      /// nodeState
      /// @param nodeState The hash map where the new states will are stored.
      /// @param expandIfNew Gives the state to set if @c index path is not 
      /// found in the map. If @c true the node is expanded. If @c false the node
      /// is collapsed.
      /// @param index Node the modify. All its children states are modified too.
      /// @param list List where all node paths the recursivity has treated are
      /// stored.
      /// @note The @c list parameter has been added to fix an endless loop bug.
      /// It @b should be removed as soon as possible if an alternative is found.
      /// @warning When calling the method for the first in the recursivity loop
      /// make sure that the list is empty to avoid incomplete read of the index
      /// and its children.
      void setNodeState(const QHash<QString, bool> & nodeState, 
                        bool expandIfNew, const QModelIndex & index, 
                        QList<QString> & list);
      
      void enableDisableOptions(const QModelIndex & index);
      
      
    }; // class TreeView
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_TreeView_h
