#include <QtXml>
#include <QDomNamedNodeMap>

#include "ClientServer/treeview/TreeItem.hh"

using namespace COOLFluiD::treeview;

TreeItem::TreeItem(QDomNode & node, TreeItem * parent)
{
  this->domNode = node;
  this->parentItem = parent;
  
  this->addChildren(node.childNodes().count());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeItem::~TreeItem()
{
  /* QHash<int,TreeItem *>::iterator it = this->childItems.begin();
   
   while(it != this->childItems.end())
   {
   delete it.value();
   it++;
   }*/
  
  if(this->parentItem != NULL)
    this->parentItem->removeChild(this);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomNode TreeItem::getDomNode() const
{
  return this->domNode;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeItem * TreeItem::getParentItem() const
{
  return this->parentItem;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeItem * TreeItem::getChild(int i)
{
  TreeItem * item = NULL;
  // if the TreeItem corresponding to this child has already been created,
  // it is returned...
  if (i>= 0 && i < this->childItems.size())
    item = childItems.at(i);
  
  // ...otherwise, if the index is valid, it is created and returned...
  if(item == NULL && i>= 0 && i < this->childItems.size())
  {
    QDomNode childNode = this->domNode.childNodes().item(i);
    item = new TreeItem(childNode, this);
    this->childItems.replace(i, item);
  }
  
  // ...if the index is not valid, return a NULL pointer
  return item;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int TreeItem::getRowNumber() 
{
  if(this->parentItem == NULL)
    return 0;
  
  return this->parentItem->getChildNumber(this);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeItem::deleteChild(int row, bool remove)
{
  if(row >= 0 && row < childItems.size())
  {
    TreeItem * item = this->childItems.at(row);
    
    if(remove)
    {
      QDomNode node = item->getDomNode();
      
      this->childItems.removeAt(row);
      node.parentNode().removeChild(node);
    }
    else
      this->childItems.replace(row, NULL); // replace by a NULL pointer
    
    delete item;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeItem::deleteChildren(int row, unsigned int count, bool remove)
{
  bool success = false;
  
  if(row >= 0 && row + count <= this->childItems.size())
  {
    QDomNodeList childNodes = this->domNode.childNodes();
    
    for(int i = 0 ; i < count ; i++)
    {
      this->deleteChild(row, remove);
      //    delete this->childItems.at(row);
      // //    this->childItems.replace(row, NULL);
    }
    
    success = true;
  }
  
  return success;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeItem::addChild(QDomNode & node)
{
  bool success = !node.isNull();
  
  if(success)
    this->childItems << new TreeItem(node, this);
  
  return success;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeItem::removeChild(TreeItem * child)
{
  return this->childItems.removeAll(child) != 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeItem::addChildren(int count)
{
  for(int i = 0 ; i < count ; i++)
    this->childItems << NULL;
  
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int TreeItem::getChildCount() const
{
  return this->childItems.count();
}

/****************************************************************************
 
 PRIVATE METHODS
 
 ****************************************************************************/

int TreeItem::getChildNumber(TreeItem * child) 
{
  return this->childItems.indexOf(child); 
}
