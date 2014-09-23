#ifndef COOLFluiD_treeview_TreeItem_h
#define COOLFluiD_treeview_TreeItem_h

/////////////////////////////////////////////////////////////////////////////

#include <QList>

#include "ClientServer/treeview/TreeViewAPI.hh"

class QDomNode;

namespace COOLFluiD
{
  namespace treeview
  {
    
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief TreeItem class represents an item of the tree model.
    
    /// @author Quentin Gasper.
    
    class TreeView_API TreeItem
    {
    private:
      
      /// @brief The node represented by this item.
      QDomNode domNode;
      
      /// @brief Children of this item.
      
      QList<TreeItem *> childItems;
      
      /// @brief A pointer to the parent item.
      
      /// This pointer may be null, if this item is the root of the tree.
      TreeItem * parentItem;
      
      int getChildNumber(TreeItem * child);
      
    public:
      /// @brief Consructor.
      
      /// @param node The node this item represents.
      /// @param parent A pointer to the parent. May be null of this item
      /// is the root of the tree.
      TreeItem(QDomNode & node, TreeItem * m_parent = NULL);
      
      /// @brief Destructor.
      
      /// Free all allocated memory.
      ~TreeItem();
      
      /// @brief Gives the child having the given row number.
      
      /// @param i Row number of the wanted child
      /// @return Returns the corresponding TreeItem, or a null pointer if the 
      /// row number is not valid.
      TreeItem * getChild(int i);
      
      /// @brief Gives the parent item of this item.
      
      /// @return Returns the parent item.
      TreeItem * getParentItem() const;
      
      /// @brief Gives the node of this item.
      
      /// @return Returns the node.
      QDomNode getDomNode() const;
      
      /// @brief Gives the row number of this item.
      
      /// @return Returns the row number.
      int getRowNumber();
      
      /// @brief Deletes a child index.
      
      /// Calling @c #getChild() after this method will recreate the index. This 
      /// is useful if a node has been replaced by another one, to avoid the
      /// m_view to crash on update (because index points to data that do not exist 
      /// anymore). If @c row does not point to a valid child, nothing is done.
      /// @param row Child row.
      void deleteChild(int row, bool remove = false);
      
      bool deleteChildren(int row, unsigned int m_count, bool remove = false);
      
      bool addChild(QDomNode & node);
      
      bool addChildren(int m_count);
      
      bool removeChild(TreeItem * child);
      
      int getChildCount() const;
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace treeview
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_treeview_TreeItem_h
