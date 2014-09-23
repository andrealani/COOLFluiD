#ifndef COOLFLuiD_treeview_TreeModel_h
#define COOLFLuiD_treeview_TreeModel_h

/////////////////////////////////////////////////////////////////////////////

#include <QAbstractItemModel>
#include <QDomDocument>
#include <QHash>
#include <QStringList>

#include "ClientServer/treeview/TreeViewAPI.hh"

class QModelIndex;
class QVariant;

namespace COOLFluiD
{
  namespace treeview
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    struct TObjectProperties;
    struct TSshInformation;
    class TreeItem;
    
    /// @brief This class provides tools to manipulate an XML tree.
    
    /// This class also provides a model (inherits 
    /// @c QAbstractItemModel class) that can be used to display the 
    /// tree in a graphical m_view using the "Model/View Programming" concept. 
    
    /// @author Quentin Gasper.
    
    class TreeView_API TreeModel : public QAbstractItemModel
    {
      Q_OBJECT
      
    private:
      /// @brief The tree
      QDomDocument m_domDocument;
      
      QStringList m_columns;
      
      /// @brief Root of the tree (used to display the tree in a m_view)
      TreeItem * m_rootItem;
      
      /// @brief Indicates wether the model is in advanced mode.
      
      /// If @c true, the model is in advanced mode. See data() for
      /// further infomation.
      bool m_advancedMode;
      
      QModelIndex currentIndex;
      
      QModelIndex currentSimulation;
      
      /// @brief Recursive method that builds a @c QStingList with 
      /// all parent nodes names of a given node. 
      
      /// The first string is the root of the tree.
      /// @param node Node from which the parents are returned
      /// @return Returns the built list
      QStringList getParentNodeNames(const QDomNode & node) const;
      
      bool isSimulationNode(const QDomNode & node) const;
      
      /// @brief Appends to a document given m_options of a node.
      
      /// @param tagName Tag name. "modOptions" for modified m_options and 
      /// "addOptions" for new m_options.
      /// @param parent Parent of the m_options (used to know the path)
      /// @param m_options Options to append.
      /// @param doc A reference to the document to which the m_options will
      /// be appended.
      /// @param keepAttrs If @c true, XML attributes are kept. 
      /// Otherwise they are removed.
      void buildModification(const QString & tagName,const QDomNode & m_parent,
                             const QDomDocument & m_options, QDomDocument & doc,
                             bool keepAttrs);
      
      /// @brief Returns the model index associated to a node.
      
      /// This is a recursive method. It calls after removing the first item 
      /// from the list until the list is empty. In each iteration, a child
      /// of the given index having the name given by the first item of the list
      /// is searched. If a such element is found, the first list item is removed
      /// the element becomes the new index. Otherwise, the recursion loop is 
      /// stopped and the method returns an invalid index.
      /// @param path Path to the wanted node. The first item is the next parent 
      /// to look for. The last item is the wanted node.
      /// @param index Current index the recursion loop is working with. The first
      /// list item should name a child of this index. When calling the method 
      /// for the first time of the recursion, root index should be passed. 
      /// @return Returns the index found, or an invalid index if the node does
      /// not exist in this tree.
      QModelIndex getIndex(QStringList & path, const QModelIndex & index) const;
      
      /// @brief Gives the index associated to the provided node
      
      /// @param node Node
      /// @return Returns the index associated to the provided node or an invalid
      /// index if the node was not found in the model.
      QModelIndex getIndex(const QDomNode & node) const;
      
      TreeItem * indexToItem(const QModelIndex & index) const;
      
      void addOption(const QString & m_name, bool basic, const QString & m_type, 
                     const QString & value, const QString & descr, 
                     QDomElement & m_parent);
      
    public:
      /// Constructor.
      
      /// @param document XML document on which this model is based.
      /// @param parent Parent of this model.
      TreeModel(QDomDocument document, QObject * m_parent = NULL);
      
      /// Destructor.
      ~TreeModel();
      
      /// @brief Implementation of @c QAbstractItemModel::data().
      
      /// Only the role @c Qt::DisplayRole is accepted. Other 
      /// roles will result to the return of an empty QVariant object 
      /// (built with the default construtor). 
      /// @param index Concerned item index.
      /// @param role Role of the returned value (only 
      /// @c Qt::DisplayRole).
      /// @return Returns an empty QVariant object if the role is not
      /// @c Qt::DisplayRole or if the @c index.isValid()
      /// returns @c false. Otherwise, returns the nodename of the
      /// the item at the specified index.
      virtual QVariant data(const QModelIndex & index, int m_role) const;
      
      /// @brief Implementation of @c QAbstractItemModel::index().
      
      /// Gives the index of the item at the given row and column under
      /// the given parent. If the parent index is not valid, the root item
      /// is taken as parent.
      /// @param row Item row from the parent.
      /// @param column Item column.
      /// @param parent Item parent.
      /// @return Returns the requested index, or a null index if
      /// <code>hasIndex(row, column, parent)</code> returns @c false.
      virtual QModelIndex index(int row, int column, 
                                const QModelIndex & m_parent = QModelIndex()) const;
      
      /// @brief Implementation of @c QAbstractItemModel::parent().
      
      /// @param child Item index of which we would like to know the parent.
      /// @return Returns the parent index of the given child or a null 
      /// index if the child is not a valid index.
      virtual QModelIndex parent(const QModelIndex &child) const;
      
      /// @brief Implementation of @c QAbstractItemModel::rowCount().
      
      /// If the parent index is not valid, the root item is taken as parent.
      /// @return Returns the row count (number of children) of a given parent.
      virtual int rowCount(const QModelIndex & m_parent = QModelIndex()) const;
      
      /// @brief Implementation of @c QAbstractItemModel::columnCount().
      /// @return Always returns 1.
      virtual int columnCount(const QModelIndex & m_parent = QModelIndex()) const;
      
      /// @brief Gives header titles.
      
      /// Overrides @c QAbstractItemModel::headerData().
      /// @param section Section number.
      /// @param orientation Header orientation.
      /// @param role Data role. Only @c Qt::DisplayRole is accepted.
      /// @return Returns the data or an empty @c QVariant on error.
      virtual QVariant headerData(int section, Qt::Orientation orientation, 
                                  int m_role = Qt::DisplayRole) const;
      
      virtual bool removeRows(int row, int m_count, const QModelIndex & index = QModelIndex());
      
      /// @brief Gives the child nodes of the given index.
      
      /// @param index The parent index.
      /// @return Returns a @c QDomNodeList containing the children,
      /// or an empty list if the provided index is not valid.
      QDomNodeList getOptions(const QModelIndex & index) const;
      
      /// @brief Builds a document that can be used as data for "modifyOption" 
      /// action.
      
      /// @param index Index from which m_options are taken.
      /// @param m_options Document cotaining modified m_options. These m_options 
      /// must be at the root of the root of the document.
      /// @param m_newOptions Document cotaining new m_options. These m_options must 
      /// be at the root of the root of the document.
      /// @return Returns the built document.
      QDomDocument modifyToDocument(const QModelIndex & index,
                                    const QDomDocument m_options,
                                    const QDomDocument m_newOptions =
                                    QDomDocument());
      
      /// @brief Gives the node associated to the given index.
      
      /// @param index Index.
      /// @return Returns the node associated, or a null node if the given
      /// index is not valid.
      QDomNode indexToNode(const QModelIndex & index) const;
      
      /// @brief Builds a node that can be used as data for "addNode" 
      /// action.
      
      /// @param index Index of the parent node.
      /// @param newNode New node name.
      /// @param doc The document the node will be added to. The presence of
      /// this parameter is due to the fact that a node can not exist if it 
      /// does not belong to a document.
      /// @return Returns the built node.
      QDomNode newChildToNode(const QModelIndex & index, 
                              const QString & newNode, QDomDocument & doc);
      
      /// @brief Builds a node that can be used as data for "renameNode" 
      /// action.
      
      /// @param index Index of the node to rename.
      /// @param newName New name of the node.
      /// @param doc The document the node will be added to. The presence of
      /// this parameter is due to the fact that a node can not exist if it 
      /// does not belong to a document.
      /// @return Returns the built node.
      QDomNode renameToNode(const QModelIndex & index,
                            const QString & newName, QDomDocument & doc);
      
      /// @brief Sets or resets the advanced mode.
      
      /// @param advanced Advanced mode state.
      void setAdvancedMode(bool advanced);
      
      /// @brief Gives the advanced mode state
      
      /// @return Returns the advanced mode state.
      bool isAdvancedMode() const;
      
      /// @brief Give the properties (XML attributes) of a given item.
      
      /// @param index Item index.
      /// @param ok Reference to a bool variable. After this method returns, 
      /// the bool value is @c false if there was an error, otherwise
      /// value is @c true.
      /// @return Returns the properties of the item. Properties values are
      /// undefined if there was an error (@c ok = @c false ).
      TObjectProperties getProperties(const QModelIndex & index, 
                                      bool & ok) const;
      
      /// @brief Gives a deep copy of the model document.
      
      /// @return Returns a deep copy of the model document.
      QDomDocument getDocument() const;
      
      /// @brief Builds a string representing the path from the root node to
      /// the node pointed by @c index.
      
      /// Path is structured like a UNIX-like path. For example : /path/to/node
      /// @param index Index for which the path has to be built.
      /// @return Returns the built path, or a null string (built with @c QString 
      /// default constructor) if @c index.isValid() returns @c false.
      QString getNodePath(const QModelIndex & index) const;
      
      /// @brief Gives the index of the node found at the given path.
      
      /// Path is structured like a UNIX-like path. For example : /path/to/node
      /// The leading @e slash is optional.
      /// @param path Path of the wanted node.
      /// @return Returns the index to that node or a non-valid index (built with
      /// @c QModelIndex if path does not exist).
      QModelIndex getIndex(const QString & path) const;
      
      bool isSimulationNode(const QModelIndex & index) const;
      
      /// @brief Creates a new simulation node
      
      /// @param simName Name of the new simulation.
      /// @return Returns @c true if the simulation was successfully created, 
      /// @c false if the creation failed (i.e. another simulation with the same 
      /// name already exists or the name is empty).
      bool createSimulation(const QString & simName);
      
      /// @brief Removes a simulation from the workspace.
      
      /// If @c index is not a simulation node, nothing is done.
      /// @param index Simulation index
      void removeSimulation(const QModelIndex & index);
      
      /// @brief Gives of a simulation with a given name.
      
      /// @param simName Simulation name.
      /// @return Returns the simulation index, or a invalid index if the 
      /// simulation was not found or if the name was empty.
      QModelIndex getSimulationIndex(const QString & simName) const;
      
      /// @brief Checks whether a simulation with a given name exists.
      
      /// @param simName Simulation name.
      /// @return Returns @c true if the simulation was found; otherwise, returns
      /// @c false.
      bool simulationExists(const QString & simName) const;
      
      /// @brief Sets simulation connection information.
      
      /// If @c index is not valid or does not point to a simulation node, 
      /// nothing is done.
      /// @param m_options Connection m_options
      /// @param index Simulation index
      void setSimConnectionInfos(const QDomDocument & m_options, 
                                 const QModelIndex & index);
      
      /// @brief Sets a simulation tree.
      
      /// If @c index is not valid or does not point to a simulation node, 
      /// nothing is done. If a tree is already present for this simulation, it
      /// is replaced.
      /// @param tree Simulation tree.
      /// @param index Simulation index.
      void setSimulationTree(const QDomDocument & tree, const QModelIndex & index);
      
      /// @brief Gives simulation connection information.
      
      /// If @c index is not valid or does not point to a simulation node, 
      /// nothing is done and @c info is not modified.
      /// @param index Simulation index.
      /// @param index Reference to a variable where information will be written.
      /// @return Returns @c true if information were successfully written into
      /// @c info ; otherwise, returns @c false.
      bool getConnectionInfos(const QModelIndex & index, 
                              TSshInformation & info) const;
      
      /// @brief Gives information about the workers.
      
      /// If @c index is not valid or does not point to a simulation node, 
      /// nothing is done and @c nbProcs and @c hosts are not modified.
      /// @param index Simulation index.
      /// @param nbProcs Reference to a variable where the number of workers
      /// will be written.
      /// @param hosts Reference to a variable where selected host machine list
      /// will be written. The variable is cleared before first write.
      /// @return Returns @c true if information were successfully written into
      /// @c nbProcs and @c hosts ; otherwise, returns @c false.
      bool getWorkersInfo(const QModelIndex & index, int & nbProcs, 
                          QStringList & m_hosts) const;
      
      /// @brief Gives parent simulation index
      
      /// This method can be used to retrieve the simulation index to which the
      /// node pointed by @c index belongs.
      /// @param index Node index
      /// @return Returns the simulation index, or an invalid index if the 
      /// simulation was not found or if @c index is not valid.
      QModelIndex getParentSimIndex(const QModelIndex & index) const;
      
      /// @brief Sets whether a simulation is connected to its server or not.
      
      /// If @c index is not valid or does not point to a simulation node, 
      /// nothing is done.
      /// @param index Simulation index.
      /// @param status If @c true, the simulation is set to "connected"; 
      /// otherwise, it is set to "not connected".
      void setSimConnectedStatus(const QModelIndex & index, bool status);
      
      /// @brief Checks wether a simulation is connected or not.
      
      /// If @c index is not valid or does not point to a simulation node, 
      /// nothing is done.
      /// @param index Simulation index
      /// @return Returns @c true if the simulation is connected; otherwise, 
      /// returns @c false.
      bool isSimulationConnected(const QModelIndex & index) const;
      
      /// @brief Sets simulation active status
      
      /// If @c index is not valid or does not point to a simulation node, 
      /// nothing is done.
      /// @param index Simulation index.
      /// @param active @c true if the simulation is active; otherwise, @c false.
      void setSimActiveState(const QModelIndex & index, bool m_active);
      
      void setCurrentSimulation(const QModelIndex & index);
      
      void setCurrentIndex(const QModelIndex & index);
      
      QModelIndex getCurrentSimulation() const;
      
      QModelIndex getCurrentIndex() const;
      
      void setSimulationHostList(const QModelIndex & index, const QDomDocument & m_hosts);
      
      QString getCurrentIndexPath() const;
      
      bool isSimulationActive(const QModelIndex & index) const;
      
      QString getCurrentPath() const;
      
      QString getSimulationName(const QModelIndex & index) const;
      
      void setSimReadOnly(const QModelIndex & index, bool m_readOnly);
      
      bool isSimReadOnly(const QModelIndex & index) const;
      
      bool isCurrentSimIndex(const QModelIndex & index) const;
      
      bool isCurrentIndex(const QModelIndex & index) const;
      
      bool areEqual(const QModelIndex & left, const QModelIndex & right) const;
      
    signals:
      
      void currentIndexChanged(const QModelIndex & index);
      
      void currentSimulationChanged(const QModelIndex & index);
      
      void advancedModeChanged(bool advanced);
      
      void readOnlyModeChanged(const QModelIndex & index, bool m_readOnly);
      
      void simulationRemoved(const QModelIndex & index);
      
    }; // class TreeModel
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace treeview
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFLuiD_treeview_TreeModel_h
