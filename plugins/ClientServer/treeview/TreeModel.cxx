#include <iostream>
#include <iterator>

#include <QtCore>
#include <QtGui>
#include <QtXml>

#include <boost/filesystem/path.hpp>

#include "Config/ConfigArgs.hh"
#include "Config/ConfigFileReader.hh"

#include "ClientServer/treeview/TreeItem.hh"
#include "ClientServer/treeview/TSshInformation.hh"
#include "ClientServer/treeview/TObjectProperties.hh"
#include "ClientServer/treeview/TreeModel.hh"

using namespace COOLFluiD::Config;
using namespace COOLFluiD::treeview;

TreeModel::TreeModel(QDomDocument document, QObject * parent)
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
    QDomElement element = m_domDocument.createElement("Workspace");
    
    element.setAttribute("tree", "object");
    element.setAttribute("mode", "basic");
    element.setAttribute("dynamic", "static");
    
    m_domDocument.appendChild(element);
  }
  
  QDomNode node = m_domDocument.firstChild();
  
  m_rootItem = new TreeItem(node, NULL);
  
  m_advancedMode = false;
  
  m_columns << "Simulations" << "Connected" << "Active";
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeModel::~TreeModel()
{
  delete m_rootItem;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int TreeModel::columnCount(const QModelIndex & parent) const
{
  return m_columns.size();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QVariant TreeModel::headerData(int section, Qt::Orientation orientation, 
                               int role) const
{
  if(role == Qt::DisplayRole && orientation == Qt::Horizontal && section >= 0 
     && section < m_columns.size())
    return m_columns.at(section);
  
  return QVariant();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::removeRows(int row, int count, const QModelIndex & index)
{
  TreeItem * parentItem;
  bool success = false;
  
  if(index.isValid())
    parentItem = this->indexToItem(index);
  else
    parentItem = m_rootItem;
  
  this->beginRemoveRows(index, row, row + count - 1);
  
  success = parentItem->deleteChildren(row, count, true);
  this->endRemoveRows();
  
  return success;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QVariant TreeModel::data(const QModelIndex & index, int role) const
{
  TreeItem *item;
  QDomNode node;
  
  if (!index.isValid() || role != Qt::DisplayRole)
    return QVariant();
  
  item = this->indexToItem(index);
  
  if(item == NULL)
    return QVariant();
  
  node = item->getDomNode();
  
  if(index.column() == 0)
  {
    QDomNamedNodeMap attributes = node.attributes();
    
    if(this->isSimulationNode(node))
      return attributes.namedItem("name").nodeValue();
    
    if(attributes.namedItem("tree").nodeValue() == "object")
    {
      if(!m_advancedMode && attributes.namedItem("mode").nodeValue() == 
         "advanced")
        return QVariant();
      
      if(!this->isSimulationNode(node))
      {
        QModelIndex parentIndex = this->getParentSimIndex(index);
        
        if(!this->isSimulationConnected(parentIndex))
          return QVariant();
      }
      
      return QString("%1 [%2]").arg(node.nodeName()).arg(attributes.namedItem("type").nodeValue());
    }
    return QVariant();
  }
  else if(this->isSimulationNode(node))
  {
    QDomNamedNodeMap attrs = node.attributes();
    
    if(index.column() == 1)
    {
      if(attrs.namedItem("connected").nodeValue() == "true")
        return "Yes";
      else
        return "No";
    }
    else if(index.column() == 2)
    {
      if(attrs.namedItem("active").nodeValue() == "true")
        return "Yes";
      else
        return "No";
    }
  }
  else
    return QVariant();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::index(int row, int column,
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
    parentItem = this->indexToItem(parent);
  
  childItem = parentItem->getChild(row);
  
  if(childItem != NULL)
    index = createIndex(row, column, childItem);
  else
    index = QModelIndex();
  
  return index;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::parent(const QModelIndex & child) const
{
  if(!child.isValid())
    return QModelIndex();
  
  TreeItem * childItem = this->indexToItem(child);
  TreeItem * parentItem = childItem->getParentItem();
  
  if (parentItem == NULL || parentItem == m_rootItem)
    return QModelIndex();
  
  return createIndex(parentItem->getRowNumber(), 0, parentItem);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int TreeModel::rowCount(const QModelIndex &parent) const
{
  TreeItem *parentItem;
  
  if (parent.column() > 0)
    return 0;
  
  if (!parent.isValid())
    parentItem = m_rootItem;
  else
    parentItem = this->indexToItem(parent);
  
  return parentItem->getDomNode().childNodes().count();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomNodeList TreeModel::getOptions(const QModelIndex & index) const
{
  TreeItem * item;
  
  if(!index.isValid())
    return QDomNodeList();
  
  item = this->indexToItem(index);
  
  if(item == NULL)
    return QDomNodeList();
  
  return item->getDomNode().childNodes();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomDocument TreeModel::modifyToDocument(const QModelIndex & index,
                                         const QDomDocument options,
                                         const QDomDocument newOptions)
{
  QDomDocument doc;
  TreeItem * item;
  
  if(!index.isValid())
    return QDomDocument();
  
  item = this->indexToItem(index);
  
  if(item == NULL)
    return QDomDocument();
  
  this->buildModification("modOptions", item->getDomNode(), options, doc, 
                          false);
  
  this->buildModification("addOptions", item->getDomNode(), newOptions, doc, 
                          true);
  
  return doc;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::buildModification(const QString & tagName,
                                  const QDomNode & parent,
                                  const QDomDocument & options,
                                  QDomDocument & doc, bool keepAttrs)
{
  QDomElement node;
  QStringList parents = this->getParentNodeNames(parent);
  QString parentsString = QString("/") + parents.join("/");
  QDomNodeList childNodes = options.childNodes();
  
  if(!childNodes.isEmpty())
  {
    node = doc.createElement(tagName);
    node.setAttribute("path", parentsString);
    
    for(int i = 0 ; i < childNodes.count() ; i++)
    {
      QDomElement option = doc.importNode(childNodes.item(i), true).toElement();
      
      if(option.isNull())
        continue;
      
      if(!keepAttrs)
      {
        QDomNamedNodeMap attributes = option.attributes();
        
        while(attributes.count() > 0)
          option.removeAttribute(attributes.item(0).nodeName());
      }
      
      node.appendChild(option);
    }
    
    doc.appendChild(node);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TObjectProperties TreeModel::getProperties(const QModelIndex & index,
                                           bool & ok) const
{
  TObjectProperties properties;
  TreeItem * item;
  QDomNamedNodeMap attributes;
  
  if(!index.isValid())
  {
    ok = false;
    return TObjectProperties();
  }
  
  item = this->indexToItem(index);
  
  if(item == NULL)
  {
    ok = false;
    return TObjectProperties();
  }
  
  attributes = item->getDomNode().attributes();
  
  properties.m_type = attributes.namedItem("type").nodeValue();
  properties.absType = attributes.namedItem("abstype").nodeValue();
  properties.dynamic = attributes.namedItem("dynamic").nodeValue() == "true";
  properties.basic = attributes.namedItem("mode").nodeValue() == "basic";
  
  return properties;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomNode TreeModel::newChildToNode(const QModelIndex & index,
                                   const QString & newNode,
                                   QDomDocument & doc)
{
  //QDomDocument doc;
  QDomElement lastElement;
  QDomElement elt;
  
  TreeItem * item;
  QDomNode node;
  QDomNode indexNode;
  
  if(!index.isValid())
    return QDomDocument();
  
  if(newNode.isNull() || newNode.isEmpty())
    return QDomDocument();
  
  item = this->indexToItem(index);
  indexNode = item->getDomNode();
  
  if(item == NULL)
    return QDomDocument();
  
  QStringList parents = this->getParentNodeNames(indexNode.parentNode());
  
  if(parents.count() > 0)
  {
    lastElement = doc.createElement(parents.at(0));
    doc.appendChild(lastElement);
    
    for(int i = 1 ; i < parents.count() ; i++)
    {
      QDomElement element = doc.createElement(parents.at(i));
      lastElement.appendChild(element);
      lastElement = element;
    }
    
    elt = doc.createElement(indexNode.nodeName());
    node = lastElement.appendChild(elt);
    
    QDomElement elem = doc.createElement(newNode);
    node.appendChild(elem);
    return elem;
  }
  
  return QDomElement();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomNode TreeModel::renameToNode(const QModelIndex & index,
                                 const QString & newName,
                                 QDomDocument & doc)
{
  QDomElement lastElement;
  QDomElement elt;
  
  TreeItem * item;
  QDomNode node;
  QDomNode indexNode;
  
  if(!index.isValid())
    return QDomDocument();
  
  if(newName.isNull() || newName.isEmpty())
    return QDomDocument();
  
  item = this->indexToItem(index);
  indexNode = item->getDomNode();
  
  if(item == NULL)
    return QDomDocument();
  
  QStringList parents = this->getParentNodeNames(indexNode.parentNode());
  
  if(parents.count() > 0)
  {
    lastElement = doc.createElement(parents.at(0));
    doc.appendChild(lastElement);
    
    for(int i = 1 ; i < parents.count() ; i++)
    {
      QDomElement element = doc.createElement(parents.at(i));
      lastElement.appendChild(element);
      lastElement = element;
    }
    
    elt = doc.createElement(indexNode.nodeName());
    elt.setAttribute("newName", newName);
    node = lastElement.appendChild(elt);
    return elt;
  }
  
  return QDomNode();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomNode TreeModel::indexToNode(const QModelIndex & index) const
{
  TreeItem * item;
  
  if(!index.isValid())
    return QDomNode();
  
  item = this->indexToItem(index);
  
  if(item == NULL)
    return QDomNode();
  
  return item->getDomNode();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setAdvancedMode(bool advanced)
{
  m_advancedMode = advanced;
  emit advancedModeChanged(advanced);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isAdvancedMode() const
{
  return m_advancedMode;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomDocument TreeModel::getDocument() const
{
  QDomDocument doc = m_domDocument.cloneNode(true).toDocument();
  QDomNodeList childNodes = doc.firstChild().childNodes();
  
  
  for(int i = 0 ; i < childNodes.size() ; i++)
  {
    QDomNode node = childNodes.at(i);
    QDomElement element = node.toElement();
    QDomNode hostListNode = node.firstChildElement("workersHosts");
    
    element.setAttribute("connected", "false");
    element.setAttribute("active", "false");
    element.setAttribute("readOnly", "false");
    
    if(!hostListNode.firstChild().isNull())
      hostListNode.removeChild(hostListNode.firstChild());
  }
  
  return doc;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString TreeModel::getNodePath(const QModelIndex & index) const
{
  TreeItem * item;
  
  if(!index.isValid())
    return QString();
  
  item = this->indexToItem(index);
  
  if(item == NULL || item->getDomNode() == m_domDocument.firstChild())
    return QString();
  
  return this->getParentNodeNames(item->getDomNode()).join("/").prepend("/");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString TreeModel::getCurrentPath() const
{
  return this->getNodePath(this->currentIndex);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString TreeModel::getSimulationName(const QModelIndex & index) const
{
  if(this->isSimulationNode(index))
  {
    TreeItem * item = this->indexToItem(index);
    QDomNamedNodeMap attrs = item->getDomNode().attributes();
    
    return attrs.namedItem("name").nodeValue();
  }
  
  return QString();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setSimReadOnly(const QModelIndex & index, bool readOnly)
{
  if(this->isSimulationNode(index))
  {
    TreeItem * item = this->indexToItem(index);
    QDomElement simNode = item->getDomNode().toElement();
    
    simNode.setAttribute("readOnly", readOnly ? "true" : "false");
    
    emit readOnlyModeChanged(index, readOnly);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isSimReadOnly(const QModelIndex & index) const
{
  bool readOnly = false;
  
  if(this->isSimulationNode(index))
  {
    TreeItem * item = this->indexToItem(index);
    QDomNamedNodeMap attrs = item->getDomNode().attributes();
    
    readOnly = attrs.namedItem("readOnly").nodeValue() == "true";
  }
  
  return readOnly;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::getIndex(const QString & path) const
{
  QStringList list;
  
  if(path.startsWith('/'))
    list = path.right(path.length() - 1).split("/", QString::SkipEmptyParts);
  else
    list = path.split("/", QString::SkipEmptyParts);
  
  if(list.isEmpty() || !this->simulationExists(list.at(0)))
    return QModelIndex();
  
  list.removeFirst();
  
  return this->getIndex(list, this->index(0, 0)); 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isSimulationNode(const QModelIndex & index) const
{
  TreeItem * item;
  
  if(!index.isValid())
    return false;
  
  item = this->indexToItem(index);
  
  if(item == NULL)
    return false;
  
  return this->isSimulationNode(item->getDomNode());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::createSimulation(const QString & simName) 
{ 
  QDomElement simNode = m_domDocument.createElement("SavedSimulation");
  QDomElement option;
  QString username;
  QRegExp regex("^USER=");
  QStringList environment = QProcess::systemEnvironment().filter(regex);
  
  if(simName.isEmpty() || this->simulationExists(simName))
    return false;
  
  if(environment.size() == 1)
    username = environment.at(0);
  
  simNode.setAttribute("name", simName);
  simNode.setAttribute("tree", "simulation");
  simNode.setAttribute("mode", "basic");
  simNode.setAttribute("dynamic", "false");
  simNode.setAttribute("connected", "false");
  simNode.setAttribute("active", "false");
  simNode.setAttribute("readOnly", "false");
  
  this->addOption("hostname", true, "std::string", "localhost", 
                  "Host name on which the server will run", simNode);
  
  this->addOption("workersCount", true, "unsigned int", "1", 
                  "Number of workers to spawn", simNode);
  
  this->addOption("workersHosts", true, "hostList", "", 
                  "Computers on which workers will be spawned", simNode);
  
  this->addOption("port", false, "unsigned int", "62784", 
                  "Port used for the connection", simNode);
  
  this->addOption("launchServerInstance", false, "bool", "false", 
                  "Check this to launch a new server instance", simNode);
  
  this->addOption("username", false, "std::string", username.remove("USER="), 
                  "Username used to connect to the remote machine when "
                  "launching a new server instance. Ignored if "
                  "launchServerInstance is not checked", simNode);
  
  QDomNode parent = m_domDocument.firstChild();
  int row = parent.childNodes().count();
  
  this->beginInsertRows(this->getIndex(parent), row, row);
  m_domDocument.firstChild().appendChild(simNode);
  m_rootItem->addChildren(1);
  this->endInsertRows();
  
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::removeSimulation(const QModelIndex & index)
{
  if(this->isSimulationNode(index))
    this->removeRows(index.row(), 1, index.parent());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::getSimulationIndex(const QString & simName) const
{
  QDomNodeList childNodes = m_domDocument.firstChild().childNodes();
  bool found = false;
  QModelIndex index;
  
  for(int i = 0 ; i < childNodes.count() && !found ; i++)
  {
    QDomNode node = childNodes.at(i);
    found = node.attributes().namedItem("name").nodeValue() == simName;
    
    if(found)
      index = this->index(i, 0);
  }
  
  return index;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::simulationExists(const QString & simName) const
{
  return this->getSimulationIndex(simName).isValid();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setSimConnectionInfos(const QDomDocument & options,
                                      const QModelIndex & index)
{
  QDomNode node = this->indexToNode(index);
  QStringList optionList;
  QStringList::iterator it;
  
  optionList << "hostname" << "username" << "port" << "launchServerInstance" 
  << "workersCount" << "workersHosts";
  
  it = optionList.begin();
  
  while(it != optionList.end()) 
  { 
    QDomText child = node.firstChildElement(*it).firstChild().toText();
    
    if(!options.firstChildElement(*it).isNull())
      child.setNodeValue(options.firstChildElement(*it).firstChild().nodeValue());
    
    it++;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setSimulationTree(const QDomDocument & tree, 
                                  const QModelIndex & index)
{
  /// @todo restore current index if in this simulation tree
  
  if(this->isSimulationNode(index))
  {
    QDomNode simNode = this->indexToNode(index);
    QDomNode node = simNode.firstChildElement("Simulator");
    
    // if the "Simulator" node was found, we replace the old tree with the new one
    
    if(!node.isNull())
    {
      int rowNumber = this->rowCount(index) - 1; //simNode.childNodes().count() - 1;
      
      this->beginRemoveRows(index, rowNumber, rowNumber);
      this->indexToItem(index)->deleteChild(rowNumber);
      simNode.removeChild(node);
      this->endRemoveRows();
      
      //    emit layoutChanged();
    }
    //   else
    {
      int rowNumber = simNode.childNodes().count();
      this->beginInsertRows(index, rowNumber, rowNumber);
      simNode.appendChild(tree.firstChild().cloneNode(true));
      this->indexToItem(index)->addChildren(1);
      this->endInsertRows();
    }
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::getConnectionInfos(const QModelIndex & index, 
                                   TSshInformation & info) const
{
  bool ok = isSimulationNode(index);
  
  if(ok)
  {
    QDomNode node = this->indexToNode(index);
    
    if(!node.isNull())
    {
      info.m_hostname = node.firstChildElement("hostname").firstChild().nodeValue();
      info.username = node.firstChildElement("username").firstChild().nodeValue();
      info.port = node.firstChildElement("port").firstChild().nodeValue().toUInt();
      info.launchServer = node.firstChildElement("launchServerInstance").firstChild().nodeValue() == "true";  
    }
  }
  
  return ok;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::getParentSimIndex(const QModelIndex & index) const
{
  if(!index.isValid())
    return QModelIndex();
  
  QDomNode node = this->indexToNode(index);
  
  while(node.nodeName() != "SavedSimulation")
    node = node.parentNode();
  
  return this->getSimulationIndex(node.attributes().namedItem("name").nodeValue());
}

/****************************************************************************
 
 PRIVATE METHODS
 
 ****************************************************************************/

QStringList TreeModel::getParentNodeNames(const QDomNode & node) const
{
  QDomNode parentNode = node.parentNode();
  QStringList list;
  
  if(parentNode.parentNode().isNull()) // if the node has no parent
    return list;
  else
  {
    list = this->getParentNodeNames(parentNode);
    
    if(this->isSimulationNode(node))
      list << node.attributes().namedItem("name").nodeValue();
    else
      list << node.nodeName();
    
    return  list;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isSimulationNode(const QDomNode & node) const
{
  return node.attributes().namedItem("tree").nodeValue() == "simulation";
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::getIndex(QStringList & path, const QModelIndex & index) const
{
  int childNumber = 0;
  QDomNodeList childNodes;
  QDomNode node;
  
  if(path.isEmpty())
    return index;
  
  childNodes = ((TreeItem *)index.internalPointer())->getDomNode().childNodes();
  
  while(childNumber < childNodes.count() && node.nodeName() != path.at(0))
  {
    node = childNodes.at(childNumber);
    if(node.nodeName() != path.at(0))
      childNumber++;
  }
  
  if(node.nodeName() == path.at(0))
  {
    QModelIndex nodeIndex = index.child(childNumber, 0);
    path.removeFirst();
    return this->getIndex(path, nodeIndex);
  } 
  else
    return QModelIndex();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::getIndex(const QDomNode & node) const
{
  return this->getIndex(this->getParentNodeNames(node).join("/").prepend('/'));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeItem * TreeModel::indexToItem(const QModelIndex & index) const
{
  return static_cast<TreeItem *>(index.internalPointer());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::addOption(const QString & name, bool basic, const QString & type, 
                          const QString & value, const QString & descr, 
                          QDomElement & parent)
{
  QDomDocument doc = parent.ownerDocument();
  QDomElement  option = doc.createElement(name);
  
  option.setAttribute("tree", "option");
  option.setAttribute("dynamic", "false");
  option.setAttribute("mode", basic ? "basic" : "advanced");
  option.setAttribute("type", type);
  option.setAttribute("description", descr);
  
  if(!value.isEmpty())
    option.appendChild(doc.createTextNode(value));
  
  parent.appendChild(option);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setSimConnectedStatus(const QModelIndex & index, bool status)
{ 
  if(this->isSimulationNode(index))
  {
    QModelIndex colIndex = this->index(index.row(), 1, index.parent());
    QDomElement node = this->indexToNode(index).toElement();
    node.setAttribute("connected", status ? "true" : "false");
    emit dataChanged(colIndex, colIndex);
    
    if(!status) // if we set the simulation to "not connected" :
    {
      // 1. we set it to "inactive"
      this->setSimActiveState(index, false);
      
      // 2. we deactivate the read-only mode
      this->setSimReadOnly(index, false);
      
      // 3. we clear the host list
      this->setSimulationHostList(index, QDomDocument());
    }
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isSimulationConnected(const QModelIndex & index) const
{
  bool connected = false;
  
  if(this->isSimulationNode(index))
  {
    QDomElement node = this->indexToNode(index).toElement();
    connected = node.attributes().namedItem("connected").nodeValue() == "true";
  }
  
  return connected;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::getWorkersInfo(const QModelIndex & index, int & nbProcs, 
                               QStringList & hosts) const
{
  bool ok = false;
  
  if(this->isSimulationNode(index))
  {
    QDomNode nbProcsNode = this->indexToNode(index).firstChildElement("workersCount");
    QDomDocument hostsNode;
    QDomNodeList childNodes;
    
    hostsNode.setContent(this->indexToNode(index).firstChildElement("workersHosts").firstChild().nodeValue());
    childNodes = hostsNode.firstChild().childNodes();
    
    hosts.clear();
    nbProcs = nbProcsNode.firstChild().nodeValue().toInt();
    
    for(int i = 0 ; i < childNodes.count() ; i++)
    {
      QDomElement node = childNodes.at(i).toElement();
      
      if(node.attribute("selected") == "true")
        hosts << node.attribute("name");
    }
    
    ok = true;
  }
  
  return ok;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setSimActiveState(const QModelIndex & index, bool active)
{
  if(this->isSimulationNode(index))
  {
    QModelIndex colIndex = index.sibling(index.row(), 2);
    QDomElement node = this->indexToNode(index).toElement();
    node.setAttribute("active", active ? "true" : "false");
    emit dataChanged(colIndex, colIndex); 
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setCurrentSimulation(const QModelIndex & index)
{
  if(this->isSimulationNode(index))
    this->currentSimulation = index;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setCurrentIndex(const QModelIndex & index)
{
  if(index.isValid())
  {
    if(this->isSimulationNode(index))
      this->currentSimulation = index;
    
    QModelIndex parentSimIndex = this->getParentSimIndex(index);
    
    this->currentIndex = index;
    
    if(parentSimIndex != this->currentSimulation)
      this->currentSimulation = parentSimIndex;
  } 
  else
  {
    this->currentSimulation = index;
    this->currentIndex = index;
  }
  
  emit currentIndexChanged(index);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::getCurrentSimulation() const
{
  return this->currentSimulation;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex TreeModel::getCurrentIndex() const
{
  return this->currentIndex;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isSimulationActive(const QModelIndex & index) const
{
  TreeItem * item = this->indexToItem(index); 
  bool active = false;
  
  if(item != NULL)
  {
    QDomNamedNodeMap attrs = item->getDomNode().attributes();
    active = attrs.namedItem("active").nodeValue() == "true";
  }
  
  return active;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeModel::setSimulationHostList(const QModelIndex & index, 
                                      const QDomDocument & hosts)
{
  if(this->isSimulationNode(index))
  {
    QDomNode simNode = this->indexToNode(index);
    QDomElement node = simNode.firstChildElement("workersHosts");
    QModelIndex optionIndex = this->getIndex(node);
    
    if(!node.firstChild().isNull())
      node.firstChild().setNodeValue(hosts.toString());
    else
      node.appendChild(m_domDocument.createTextNode(hosts.toString()));
    
    emit dataChanged(index , index);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isCurrentSimIndex(const QModelIndex & index) const
{
  return index.isValid() && index == this->currentSimulation;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::isCurrentIndex(const QModelIndex & index) const
{
  return index.isValid() && index == this->currentIndex;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeModel::areEqual(const QModelIndex & left, const QModelIndex & right) const
{
  return left.internalPointer() == right.internalPointer() && left.row() == right.row();
}
