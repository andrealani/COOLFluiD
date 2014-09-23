#ifndef COOLFLuiD_treeview_StatusModel_h
#define COOLFLuiD_treeview_StatusModel_h

/////////////////////////////////////////////////////////////////////////////

#include <QAbstractItemModel>
#include <QDomDocument>
#include <QList>
#include <QStringList>

#include "ClientServer/treeview/TreeItem.hh"


class QModelIndex;
class QVariant;

namespace COOLFluiD
{
  namespace client
  {
    class StatusModelSubsystem;
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class provides tools to manipulate an XML tree to manage 
    /// workers status.
    
    /// This class also provides a model (inherits 
    /// @c QAbstractItemModel class) that can be used to display the 
    /// tree in a graphical m_view using the "Model/View Programming" concept. 
    
    /// @author Quentin Gasper.
    
    class StatusModel : public QAbstractItemModel
    {
      Q_OBJECT
      
    public:
      /// @brief Constructor.
      
      /// @param document XML document on which this model is based.
      /// @param parent Parent of this model.
      StatusModel(QDomDocument document, QObject * parent = NULL);
      
      /// @brief Destructor.
      ~StatusModel();
      
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
      virtual QVariant data(const QModelIndex & index, int role) const;
      
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
                                const QModelIndex & parent = QModelIndex()) const;
      
      /// @brief Implementation of @c QAbstractItemModel::parent().
      
      /// @param child Item index of which we would like to know the parent.
      /// @return Returns the parent index of the given child or a null 
      /// index if the child is not a valid index.
      virtual QModelIndex parent(const QModelIndex &child) const;
      
      /// @brief Implementation of @c QAbstractItemModel::rowCount().
      
      /// If the parent index is not valid, the root item is taken as parent.
      /// @return Returns the row count (number of children) of a given parent.
      virtual int rowCount(const QModelIndex & parent = QModelIndex()) const;
      
      /// @brief Implementation of @c QAbstractItemModel::columnCount().
      /// @return Always returns 1.
      virtual int columnCount(const QModelIndex & parent = QModelIndex()) const;
      
      /// @brief Gives header titles.
      
      /// Overrides @c QAbstractItemModel::headerData().
      /// @param section Section number.
      /// @param orientation Header orientation.
      /// @param role Data role. Only @c Qt::DisplayRole is accepted.
      /// @return Returns the data or an empty @c QVariant on error.
      virtual QVariant headerData(int section, Qt::Orientation orientation, 
                                  int role = Qt::DisplayRole) const;
      
      /// @brief Gives a deep copy of the model document.
      
      /// @return Returns a deep copy of the model document.
      QDomDocument getDocument() const;
      
      void setWorkerStatus(const QString & subSysName, int rank, 
                           const QString & status);
      
      void addSubSystem(const QString & subSysName, unsigned int size);
      
      void removeSubSystem(const QString & subSysName);
      
      void clear();
      
    signals:
      
      void subSystemAdded(const QModelIndex & index);
      
    private:
      /// @brief The tree
      QDomDocument m_domDocument;
      
      QStringList m_columns;
      
      /// @brief Root of the tree (used to display the tree in a m_view)
      COOLFluiD::treeview::TreeItem * m_rootItem;
      
      QHash<QString, QHash<int, QString> *> m_modelData;
      
      QList<StatusModelSubsystem *> m_subSystems;
      
      bool isSubSystem(const QModelIndex & index) const;
      
      
    }; // class TreeModel
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace treeview
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFLuiD_treeview_StatusModel_h
