#include <QtGui>
#include <QtXml>
#include <QtCore>
#include <stdexcept>

#include "ClientServer/client/CommitDetails.hh"
#include "ClientServer/client/OptionPanel.hh"
#include "ClientServer/client/CommitDetailsDialog.hh"
#include "ClientServer/client/MenuActionInfo.hh"
#include "ClientServer/client/UnknownTypeException.hh"
#include "ClientServer/treeview/TObjectProperties.hh"
#include "ClientServer/treeview/TreeModel.hh"
#include "ClientServer/treeview/TreeItem.hh"
#include "ClientServer/treeview/TSshInformation.hh"
#include "ClientServer/client/ConfirmCommitDialog.hh"

#include "ClientServer/client/TreeView.hh"

using namespace COOLFluiD::client;
using namespace COOLFluiD::treeview;

TreeView::TreeView(OptionPanel * optionsPanel, QMainWindow * parent)
: QTreeView(parent)
{
  MenuActionInfo config;
  
  if(optionsPanel == NULL)
    throw std::invalid_argument("Options panel is a null pointer");
  
  // instantiate class attributes
  m_treeModel = NULL;
  m_modelFilter = new QSortFilterProxyModel();
  m_mnuNewOption = new QMenu("Add an option");
  
  m_optionsPanel = optionsPanel;
  
  this->setModel(m_modelFilter);
  m_modelFilter->setDynamicSortFilter(true);
  
  QRegExp reg(QRegExp(".+", Qt::CaseInsensitive, QRegExp::RegExp));
  m_modelFilter->setFilterRegExp(reg);
  
  this->buildSimulationMenu();
  this->buildObjectMenu();
  
  this->setReadOnly(false);
  
  // when right clicking on the treeview, 
  // a "Context menu event" must be generated
  this->setContextMenuPolicy(Qt::CustomContextMenu);
  
  this->setHeaderHidden(false);
  
  this->header()->setResizeMode(QHeaderView::ResizeToContents);
  this->header()->setStretchLastSection(false);
  
  connect(this, SIGNAL(activated(const QModelIndex &)), 
          this, SLOT(nodeActivated(const QModelIndex &)));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeView::~TreeView()
{
  delete m_simulationMenu;
  delete m_objectMenu;
  delete m_mnuAbstractTypes;
  delete m_modelFilter;
  delete m_treeModel;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::setAbstractTypesList(const QStringList & abstractTypes)
{
  QStringList::const_iterator it = abstractTypes.begin();
  
  // if the menu is not enabled, it can not be modified, even by code.
  m_mnuAbstractTypes->setEnabled(true);
  
  m_abstractTypesActions.clear();
  
  // clearing the menu will destroy all its actions (the ones that are not
  // linked to another menu, toolbar, etc...)...so no additional delete is
  // needed
  m_mnuAbstractTypes->clear();
  
  // create menu m_items for each abstract type
  while(it != abstractTypes.end())
  {
    MenuActionInfo actionInfo;
    QAction * action;
    
    actionInfo.m_text = *it;
    actionInfo.m_slot = SLOT(addNode());
    actionInfo.m_menu = m_mnuAbstractTypes;
    
    action = actionInfo.buildAction(this);
    m_abstractTypesActions.append(action);
    it++;
  }
  
  // if there is at least one abstract type, the sub-menu can be enabled
  m_mnuAbstractTypes->setEnabled(!m_abstractTypesActions.isEmpty());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomNode TreeView::newChildNode(const QString & newNode, QDomDocument & doc) const
{
  QModelIndex index = m_modelFilter->mapToSource(m_treeModel->getCurrentIndex());
  return m_treeModel->newChildToNode(index, newNode, doc);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::setReadOnly(bool readOnly)
{
  m_readOnly = readOnly;
  m_mnuAbstractTypes->setEnabled(!readOnly);
  m_mnuNewOption->setEnabled(!readOnly);
  
  m_actions[ ACTION_OBJECT_DELETE ]->setEnabled(!readOnly);
  m_actions[ ACTION_OBJECT_RENAME_NODE ]->setEnabled(!readOnly);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeView::isReadOnly() const
{
  return m_readOnly;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::addNode()
{
  QAction * action = static_cast<QAction *>(sender());
  
  if(action != NULL && m_abstractTypesActions.contains(action) && 
     this->confirmChangeOptions(m_treeModel->getCurrentIndex(), true))
    emit addNode(action->text());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::mousePressEvent(QMouseEvent * event)
{ 
  QPoint mousePosition(event->x() + this->x(), event->y() + this->y());
  QModelIndex index = this->indexAt(mousePosition);
  
  QModelIndex indexInModel = m_modelFilter->mapToSource(index);
  Qt::MouseButton button = event->button();
  
  this->enableDisableOptions(m_treeModel->getParentSimIndex(indexInModel));
  
  if(event->type() == QEvent::MouseButtonDblClick && button == Qt::LeftButton
     && index.isValid())
  {
    if(this->isExpanded(index))
      this->collapse(index);
    else
      this->expand(index);
  }
  else if(button == Qt::RightButton)
  {
    if(!m_treeModel->getCurrentIndex().isValid())
      m_treeModel->setCurrentIndex(indexInModel);
    
    // if this is a simulation node, we display the simulation context menu
    if(m_treeModel->isSimulationNode(indexInModel) || !indexInModel.isValid())
      m_simulationMenu->exec(QCursor::pos());
    else // we display the objects context menu
      m_objectMenu->exec(QCursor::pos());
  }
  else if(!m_treeModel->areEqual(indexInModel, m_treeModel->getCurrentIndex()))
  {  
    if(index != m_treeModel->getCurrentIndex() && this->confirmChangeOptions(index))
      m_treeModel->setCurrentIndex(indexInModel);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::keyPressEvent(QKeyEvent * event)
{
  if(event->key() == Qt::Key_Up)
  {
    QModelIndex index = m_treeModel->getCurrentIndex();
    
    if(this->confirmChangeOptions(index, true))
    {
      QModelIndex above = this->indexAbove(m_modelFilter->mapFromSource(index));
      
      if(above.isValid())
        m_treeModel->setCurrentIndex(m_modelFilter->mapToSource(above));
    }
  }
  else if(event->key() == Qt::Key_Down)
  {
    QModelIndex index = m_treeModel->getCurrentIndex();
    
    if(this->confirmChangeOptions(index, true))
    {
      QModelIndex below = this->indexBelow(m_modelFilter->mapFromSource(index));
      
      if(below.isValid())
        m_treeModel->setCurrentIndex(m_modelFilter->mapToSource(below));
    }  
  }
  else
    QTreeView::keyPressEvent(event); 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::buildObjectMenu()
{
  MenuActionInfo actionInfo;
  
  // context menu init
  QStringList typesList = OptionTypes::getTypesList();
  QStringList::iterator it = typesList.begin();
  
  // build menu m_items 
  actionInfo.m_slot = SLOT(addOption());
  actionInfo.m_menu = m_mnuNewOption;
  
  while(it != typesList.end())
  {
    actionInfo.m_text = *it;
    actionInfo.buildAction(this);
    it++;
  }
  
  //----------------------------------------------------
  //----------------------------------------------------
  
  m_mnuAbstractTypes = new QMenu("Add a child node");
  m_objectMenu = new QMenu("Context menu");
  
  m_objectMenu->addMenu(m_simulationMenu);
  m_objectMenu->addMenu(m_mnuAbstractTypes);
  m_objectMenu->addMenu(m_mnuNewOption);
  
  //--------------------------------------------
  
  m_objectMenu->addSeparator();
  
  //--------------------------------------------
  
  actionInfo.m_menu = m_objectMenu;
  actionInfo.m_text = "Rename";
  actionInfo.m_slot = SLOT(renameNode());
  
  m_actions[ACTION_OBJECT_RENAME_NODE] = actionInfo.buildAction(this);
  
  //--------------------------------------------
  
  m_objectMenu->addSeparator();
  
  //--------------------------------------------
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_objectMenu;
  actionInfo.m_text = "Delete";
  actionInfo.m_slot = SLOT(deleteNode());
  
  m_actions[ACTION_OBJECT_DELETE] = actionInfo.buildAction(this);
  
  //--------------------------------------------
  
  m_objectMenu->addSeparator();
  
  //--------------------------------------------
  
  actionInfo.initDefaults();
  actionInfo.m_menu = m_objectMenu;
  actionInfo.m_text = "Properties";
  actionInfo.m_slot = SLOT(showProperties());
  
  m_actions[ACTION_OBJECT_DELETE] = actionInfo.buildAction(this);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::buildSimulationMenu()
{
  MenuActionInfo config;
  m_simulationMenu = new QMenu("Simulation");
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "New simulation";
  config.m_slot = SLOT(newSimulation());
  
  m_actions[ACTION_SIM_NEW_SIMULATION] = config.buildAction(this);
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Open simulation";
  config.m_slot = SLOT(openSimulation());
  
  m_actions[ACTION_SIM_OPEN_SIMULATION] = config.buildAction(this);
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Close simulation";
  config.m_slot = SLOT(endSimulation());
  
  m_actions[ACTION_SIM_END_SIMULATION] = config.buildAction(this);
  
  //--------------------------------------------
  
  m_simulationMenu->addSeparator();
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Connect";
  config.m_slot = SLOT(connectSimulation());
  
  m_actions[ACTION_SIM_CONNECT] = config.buildAction(this);
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Disconnect";
  config.m_slot = SLOT(disconnectSimulation());
  
  m_actions[ACTION_SIM_DISCONNECT] = config.buildAction(this);
  
  //--------------------------------------------
  
  m_simulationMenu->addSeparator();
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Update tree";
  config.m_slot = SLOT(updateTree());
  
  m_actions[ACTION_SIM_UPDATE_TREE] = config.buildAction(this);
  
  //--------------------------------------------
  
  m_simulationMenu->addSeparator();
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Run simulation";
  config.m_slot = SLOT(runSimulation());
  
  m_actions[ACTION_SIM_RUN_SIMULATION] = config.buildAction(this);
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Stop simulation";
  config.m_slot = SLOT(stopSimulation());
  
  m_actions[ACTION_SIM_STOP_SIMULATION] = config.buildAction(this);
  
  //--------------------------------------------
  
  m_simulationMenu->addSeparator();
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Activate simulation";
  config.m_slot = SLOT(activateSimulation());
  
  m_actions[ACTION_SIM_ACTIVATE_SIM] = config.buildAction(this);
  
  //--------------------------------------------
  
  config.initDefaults();
  config.m_menu = m_simulationMenu;
  config.m_text = "Deactivate simulation";
  config.m_slot = SLOT(deactivateSimulation());
  
  m_actions[ACTION_SIM_DEACTIVATE_SIM] = config.buildAction(this);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::deleteNode()
{
  // get the index in the filter
  QModelIndex index = this->currentIndex();
  
  // get the corresponding index in the model
  QModelIndex indexInModel = m_modelFilter->mapToSource(index);
  QDomNode node = m_treeModel->indexToNode(indexInModel);
  
  emit deleteNode(node);
  
  // if the m_options in the m_options panel belong to the deleted node, 
  // we delete them
  //  if(index == m_treeModel->getCurrentIndex())
  //   this->optionsPanel->setOptions(QDomNodeList());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::renameNode()
{
  QModelIndex index = this->currentIndex();
  QModelIndex indexInModel = m_modelFilter->mapToSource(index);
  QDomNode node = m_treeModel->indexToNode(indexInModel);
  
  if(!node.isNull())
  {
    QString name = QInputDialog::getText(this, "Rename node",
                                         "New name of the new node:",
                                         QLineEdit::Normal, node.nodeName());
    
    // remove starting and ending spaces
    name = name.trimmed();
    // replace spaces by underscores
    name = name.replace(" ", "_");
    
    if(!name.isNull() && !name.isEmpty() && name != node.nodeName())
    {
      QDomDocument doc;
      QDomNode node2 = m_treeModel->renameToNode(indexInModel, name, doc);
      emit renameNode(node2, name);
    }
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::showProperties()
{
  bool ok;
  
  // get the index in the filter
  QModelIndex index = this->currentIndex();
  
  // get the corresponding index in the model
  QModelIndex indexInModel = m_modelFilter->mapToSource(index);
  
  TObjectProperties properties = m_treeModel->getProperties(indexInModel, ok);
  
  QString message = QString("Abstract type: %1\n").arg(properties.absType);
  message += QString("Type: %1\n").arg(properties.m_type);
  message += QString("Mode: %1\n").arg(properties.basic ? "basic" : "advanced");
  message += QString("Dynamic: %1\n").arg(properties.dynamic ? "yes" : "no");
  
  QMessageBox::information(this, "Properties", message);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::addOption()
{
  // action (menu item) on which user clicked
  QAction * action = qobject_cast<QAction *>(sender()); 
  
  QModelIndex index = m_treeModel->getCurrentIndex();
  QModelIndex oldIndex = m_treeModel->getCurrentIndex();
  TOptionTypes type = OptionTypes::getTypeId(action->text());
  
  if(this->confirmChangeOptions(index))
  {
    QString name;
    
    // change m_options if another node is selected
    if(index != m_treeModel->getCurrentIndex())
    {
      QModelIndex indexInModel = m_modelFilter->mapToSource(index);
      QDomNodeList options = m_treeModel->getOptions(indexInModel);
      
      m_treeModel->getCurrentIndex() = index;
      this->setCurrentIndex(index);
    }
    
    name = QInputDialog::getText(this->parentWidget(), "New option",
                                 "Enter the name of the new option:");
    
    
    if(!name.isNull() && !name.isEmpty())
      m_optionsPanel->addOption(type, name);
    
    // restore old m_options if user cancels 
    else 
    {
      QModelIndex indexInModel = m_modelFilter->mapToSource(oldIndex);
      QDomNodeList options = m_treeModel->getOptions(indexInModel);
      
      m_treeModel->getCurrentIndex() = oldIndex;
      this->setCurrentIndex(oldIndex);
    }
  } // "if(this->confirmChangeOptions(index))"
  
  else
    this->setCurrentIndex(m_treeModel->getCurrentIndex());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool TreeView::confirmChangeOptions(const QModelIndex & index, bool okIfSameIndex)
{
  bool confirmed = true;
  QMessageBox question(this);
  
  if(!okIfSameIndex && m_treeModel->areEqual(m_treeModel->getCurrentIndex(), index))
    return confirmed;
  
  if(m_optionsPanel->isModified())
  {
    CommitDetails commitDetails;
    
    ConfirmCommitDialog dlg;
    
    m_optionsPanel->getModifiedOptions(commitDetails);
    
    ConfirmCommitDialog::CommitConfirmation answer = dlg.show(commitDetails);
    
    if(answer == ConfirmCommitDialog::COMMIT)
      m_optionsPanel->commitChanges();
    
    confirmed = answer != ConfirmCommitDialog::CANCEL;
  }
  
  return confirmed;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::changesMade(const QDomDocument & modOptions, 
                           const QDomDocument & newOptions)
{
  QDomDocument doc;
  
  // get the index in the filter
  QModelIndex index = m_treeModel->getCurrentIndex();
  
  // get the corresponding index in the model
  QModelIndex indexInModel = m_modelFilter->mapToSource(index);
  
  if(m_treeModel->isSimulationNode(indexInModel))
    m_treeModel->setSimConnectionInfos(modOptions, indexInModel);
  else
  {
    doc = m_treeModel->modifyToDocument(indexInModel, modOptions,
                                              newOptions);
    emit commitChanges(doc);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::getChildrenNodesState(const QModelIndex & index,
                                     QHash<QString, bool> & map)
{
  if(!index.isValid())
    return;
  
  QModelIndex childIndex = index.child(0, 0);
  QModelIndex indexInModel = m_modelFilter->mapToSource(index);
  QString nodePath = m_treeModel->getNodePath(indexInModel);
  
  // if the path already exists in the map, nothing to do
  if(map.contains(nodePath))
    return;
  
  map[nodePath] = this->isExpanded(index);
  
  while(childIndex.isValid())
  {
    this->getChildrenNodesState(childIndex, map);
    childIndex = this->indexBelow(childIndex);
  }
  
  this->getChildrenNodesState(this->indexBelow(index), map);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::setNodeState(const QHash<QString, bool> & nodeState, 
                            bool expandIfNew, const QModelIndex & index,
                            QList<QString> & list)
{
  if(!index.isValid())
    return;
  
  bool expand;
  QModelIndex childIndex = index.child(0, 0);
  QModelIndex indexInModel = m_modelFilter->mapToSource(index);
  QString nodePath = m_treeModel->getNodePath(indexInModel);
  
  // bug fixing (endless loop): if the path already exists in the 
  // list, nothing to do
  if(list.contains(nodePath))
    return;
  
  list.append(nodePath);
  
  if(nodeState.contains(nodePath))
    expand = nodeState[nodePath];
  else
    expand = expandIfNew;
  
  if(expand)
    this->expand(index);
  else
    this->collapse(index);
  
  while(childIndex.isValid())
  {
    this->setNodeState(nodeState, expandIfNew, childIndex, list);
    childIndex = this->indexBelow(childIndex);
  }
  
  this->setNodeState(nodeState, expandIfNew, this->indexBelow(index), list);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::newSimulation()
{
  QString simName;
  bool stop = false;
  
  do
  {
    simName = QInputDialog::getText(NULL, "New simulation", "Please "
                                    "enter the name of the new simulation.");
    
    stop = simName.isEmpty();
    
    if(!stop)
    {
      stop = m_treeModel->createSimulation(simName);
      
      if(!stop)
        QMessageBox::critical(NULL, "Error", "A simulation named \"" + simName + 
                              "\" already exists.");
    }
  } while(!stop);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::openSimulation()
{
  emit openSimulation(m_treeModel->getCurrentSimulation());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::endSimulation()
{
  QMessageBox msgBox(this->parentWidget());
  QModelIndex index = m_treeModel->getCurrentSimulation();
  QString text = "This will remove the simulation \"%1\".";
  QString infoText = "Are you sure to continue?";
  
  msgBox.setWindowTitle("Close simulation");
  msgBox.setIcon(QMessageBox::Warning);
  msgBox.setText(text.arg(m_treeModel->getSimulationName(index)));
  msgBox.setInformativeText(infoText);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  
  if(msgBox.exec() == QMessageBox::Yes)
    m_treeModel->removeSimulation(index); 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::connectSimulation()
{ 
  QModelIndex index = m_treeModel->getCurrentIndex();
  
  if(this->confirmChangeOptions(index, false))
  {
    TSshInformation info;
    m_treeModel->getConnectionInfos(index, info);
    emit connectSimulation(m_treeModel->getCurrentSimulation(), info);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::disconnectSimulation()
{
  QModelIndex index = m_treeModel->getCurrentIndex();
  bool shutdown;
  bool disc;
  
  
  if(this->confirmChangeOptions(index, false))
  {
    QMessageBox discBox(this);
    QPushButton * btDisc = NULL;
    QPushButton * btCancel = NULL;
    QPushButton * btShutServer = NULL;
    
    btDisc = discBox.addButton("Disconnect", QMessageBox::NoRole);
    btCancel = discBox.addButton(QMessageBox::Cancel);
    btShutServer = discBox.addButton("Shutdown server", QMessageBox::YesRole);
    
    discBox.setIcon(QMessageBox::Question);
    discBox.setWindowTitle("Confirmation");
    discBox.setText("You are about to disconnect from the server.");
    discBox.setInformativeText("What do you want to do?");
    
    // show the message box
    discBox.exec();
    
    disc = discBox.clickedButton() != btCancel;
    shutdown = discBox.clickedButton() == btShutServer;
    
    delete btDisc;
    delete btCancel;
    delete btShutServer;
    
    if(disc)
      emit disconnectSimulation(m_treeModel->getCurrentSimulation(), shutdown);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::updateTree()
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::runSimulation()
{
  QModelIndex index = m_treeModel->getCurrentIndex();
  
  if(this->confirmChangeOptions(index, false))
    emit runSimulation(m_treeModel->getCurrentSimulation());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::stopSimulation()
{
  emit stopSimulation(m_treeModel->getCurrentSimulation());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::activateSimulation()
{
  QModelIndex index = m_treeModel->getCurrentSimulation();
  
  if(this->confirmChangeOptions(index, false))
    emit activateSimulation(m_treeModel->getCurrentSimulation());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::deactivateSimulation()
{
  emit deactivateSimulation(m_treeModel->getCurrentSimulation());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::currentIndexChanged(const QModelIndex & index)
{
  QItemSelectionModel::SelectionFlags flags = QItemSelectionModel::Select | QItemSelectionModel::Rows;
  QModelIndex indexInFilter = m_modelFilter->mapFromSource(index);
  
  this->selectionModel()->clearSelection();
  this->selectionModel()->select(indexInFilter, flags);
  this->selectionModel()->setCurrentIndex(indexInFilter, flags);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QAction * TreeView::addSimToMenu(QMenu * menu)
{
  QAction * action = NULL;
  
  if(menu != NULL)
    action = menu->addMenu(m_simulationMenu);
  
  return action;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QAction * TreeView::addSimToMenuBar(QMenuBar * menuBar)
{
  QAction * action = NULL;
  
  if(menuBar != NULL)
    action = menuBar->addMenu(m_simulationMenu);
  
  return action;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::enableDisableOptions(const QModelIndex & index)
{
  bool connected = m_treeModel->isSimulationConnected(index);
  bool active = m_treeModel->isSimulationActive(index);
  
  m_actions[ ACTION_SIM_OPEN_SIMULATION ]->setEnabled(m_treeModel->isSimulationNode(index));
  m_actions[ ACTION_SIM_END_SIMULATION ]->setEnabled(index.isValid() && !connected);
  m_actions[ ACTION_SIM_RUN_SIMULATION ]->setEnabled(connected);
  m_actions[ ACTION_SIM_STOP_SIMULATION ]->setEnabled(connected);
  
  m_actions[ ACTION_SIM_CONNECT ]->setEnabled(index.isValid() && !connected);
  m_actions[ ACTION_SIM_DISCONNECT ]->setEnabled(connected);
  m_actions[ ACTION_SIM_UPDATE_TREE ]->setEnabled(connected);
  
  m_actions[ ACTION_SIM_ACTIVATE_SIM ]->setEnabled(connected && !active);
  m_actions[ ACTION_SIM_DEACTIVATE_SIM ]->setEnabled(connected && active);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::setTreeModel(TreeModel * treeModel)
{
  if(m_treeModel != NULL)
    m_treeModel->disconnect();
  
  m_treeModel = treeModel;
  
  m_modelFilter->setSourceModel(m_treeModel);
  
  if(m_treeModel != NULL)
  {
    connect(m_treeModel, SIGNAL(currentIndexChanged(const QModelIndex &)), 
            this, SLOT(currentIndexChanged(const QModelIndex &)));
    
    // enable/disable menu m_options
    this->enableDisableOptions(m_treeModel->getCurrentIndex());
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TreeModel * TreeView::getTreeModel() const
{
  return m_treeModel;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TreeView::nodeActivated(const QModelIndex & index)
{
  
}
