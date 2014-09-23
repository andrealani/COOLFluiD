#include <iostream>
#include <iterator>

#include <QtGui>
#include <QtXml>
#include <QtCore>

#include <boost/filesystem/path.hpp>

#include "Common/PE.hh"

#include "Config/ConfigArgs.hh"
#include "Config/ConfigFileReader.hh"

#include "ClientServer/client/GlobalLog.hh"
#include "ClientServer/client/StatusModel.hh"

using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::treeview;
using namespace COOLFluiD::client;

StatusModel::StatusModel(QDomDocument document, QObject * parent)
: QAbstractItemModel(parent)
{
  QDomNodeList nodeList;
  m_domDocument = document;
  
  nodeList = m_domDocument.childNodes();
  
  // if the first node is the xml tag (<?xml...), it's removed (the
  // second node becomes the first one) : there's no need to show it.
  if(nodeList.item(0).nodeType() == QDomNode::ProcessingInstructionNode)
    m_domDocument.replaceChild(nodeList.item(1), nodeList.item(0));
  
  if(m_domDocument.childNodes().isEmpty())
  {
    QDomElement root = m_domDocument.createElement("ROOT");
    
    root.appendChild(m_domDocument.createElement("Workers"));
    m_domDocument.appendChild(root);
  }
  
  QDomNode node = m_domDocument.firstChild();
  m_rootItem = new TreeItem(node, NULL);
  
  m_columns << "Workers" << "Status";
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

StatusModel::~StatusModel()
{
  delete m_rootItem;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int StatusModel::columnCount(const QModelIndex & parent) const
{
  return m_columns.size();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QVariant StatusModel::headerData(int section, Qt::Orientation orientation, 
                                 int role) const
{
  if(role == Qt::DisplayRole && orientation == Qt::Horizontal && section >= 0 
     && section < m_columns.size())
    return m_columns.at(section);
  
  return QVariant();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomDocument StatusModel::getDocument() const
{
  return m_domDocument.cloneNode(true).toDocument();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void StatusModel::setWorkerStatus(const QString & subSysName, int rank,
                                  const QString & status)
{
  QDomNodeList childNodes = m_domDocument.firstChild().firstChild().childNodes();
  QDomNode node;
  int childCount = childNodes.count();
  QDomElement workerNode; 
  int rowNumber = 0;
  QModelIndex parentIndex;
  QModelIndex workerIndex;
  
  for(int i = 0 ; i < childCount && node.isNull() ; i++)
  {
    QDomNode tmpNode = childNodes.at(i);
    
    if(tmpNode.nodeName() == subSysName)
      node = tmpNode;
    else
      rowNumber++;
  }
  
  if(node.isNull() || status.isEmpty())
    return;
  
  parentIndex = this->index(rowNumber, 0, this->index(0, 0, QModelIndex()));
  childNodes = node.childNodes();
  workerNode = childNodes.at(rank).toElement();
  
  if(workerNode.isNull())
    return;
  
  workerNode.setAttribute("status", status);
  workerIndex = this->index(rank, 1, parentIndex);
  
  emit dataChanged(workerIndex, workerIndex);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void StatusModel::addSubSystem(const QString & subSysName, unsigned int size)
{
  QDomElement node = m_domDocument.firstChild().firstChild().namedItem(subSysName).toElement();
  
  if(!node.isNull() || size == 0)
    return;
  else
  {
    int rowNumber = m_domDocument.firstChild().firstChild().childNodes().count();
    QModelIndex index = this->index(0, 0, QModelIndex());
    
    node = m_domDocument.createElement(subSysName);
    
    for(unsigned int i = 0 ; i < size ; i++)
    {
      QDomElement workerNode = m_domDocument.createElement("worker");
      
      workerNode.setAttribute("rank", i);
      workerNode.setAttribute("status", WorkerStatus::Convert::to_str(WorkerStatus::NOT_RUNNING).c_str());
      node.appendChild(workerNode);
    }
    
    this->beginInsertRows(index, rowNumber, rowNumber);
    ((TreeItem*)index.internalPointer())->addChildren(1);
    m_domDocument.firstChild().firstChild().appendChild(node);
    this->endInsertRows();
    
    index = this->index(rowNumber, 0, index);
    emit subSystemAdded(index);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void StatusModel::removeSubSystem(const QString & subSysName)
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void StatusModel::clear()
{
  QModelIndex index = this->index(0, 0, QModelIndex());
  TreeItem * item = static_cast<TreeItem*>(index.internalPointer());
  
  if(index.isValid() && item != NULL)
  {
    QDomNode node = item->getDomNode().parentNode();
    QDomNodeList childNodes = node.childNodes();
    int rowCount = this->rowCount(QModelIndex());
    
    this->beginRemoveRows(QModelIndex(), 0, rowCount - 1);
    
    for(int i = 0 ; i < rowCount ; i++)
      node.removeChild(childNodes.at(i));
    
    this->endRemoveRows();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QVariant StatusModel::data(const QModelIndex & index, int role) const
{
  TreeItem *item;
  QDomNode node;
  
  if (!index.isValid() || role != Qt::DisplayRole)
    return QVariant();
  
  item = static_cast<TreeItem*>(index.internalPointer());
  
  if(item == NULL)
    return QVariant();
  
  node = item->getDomNode();
  
  QDomNamedNodeMap attributes = node.attributes();
  
  if(index.column() == 0)
  {
    if(attributes.isEmpty())
      return node.nodeName();
    else
      return QString("Worker[%1]").arg(attributes.namedItem("rank").nodeValue());
  }
  else if(index.column() == 1)
  {
    if(attributes.isEmpty())
      return QVariant();
    else
      return attributes.namedItem("status").nodeValue();
  }
  else
    return QVariant();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex StatusModel::index(int row, int column,
                               const QModelIndex &parent) const
{
  TreeItem * childItem;
  TreeItem * parentItem;
  QModelIndex index;
  
  if(!this->hasIndex(row, column, parent))
    return QModelIndex();
  
  if(!parent.isValid())
    parentItem = m_rootItem;
  else
    parentItem = static_cast<TreeItem*>(parent.internalPointer());
  
  childItem = parentItem->getChild(row);
  
  if(childItem != NULL)
    index = createIndex(row, column, childItem);
  else
    index = QModelIndex();
  
  return index;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex StatusModel::parent(const QModelIndex & child) const
{
  if(!child.isValid())
    return QModelIndex();
  
  TreeItem * childItem = static_cast<TreeItem*>(child.internalPointer());
  TreeItem * parentItem = childItem->getParentItem();
  
  if (parentItem == NULL || parentItem == m_rootItem)
    return QModelIndex();
  
  return createIndex(parentItem->getRowNumber(), 0, parentItem);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int StatusModel::rowCount(const QModelIndex &parent) const
{
  TreeItem *parentItem;
  
  if (parent.column() > 0)
    return 0;
  
  if (!parent.isValid())
    parentItem = m_rootItem;
  else
    parentItem = static_cast<TreeItem*>(parent.internalPointer());
  
  return parentItem->getDomNode().childNodes().count();
}
