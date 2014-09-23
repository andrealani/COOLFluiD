#include <QtCore>
#include <QtGui>

#include "ClientServer/client/RemoteOpenFile.hh"
#include "ClientServer/client/FilesPanel.hh"

using namespace COOLFluiD::client;

FilesPanel::FilesPanel(bool includeFiles, const QStringList & extensions, 
                       bool includeNoExtension, QWidget * parent) 
: QWidget(parent)
{ 
  //  this->openFileDialog = new RemoteOpenFile();
  
  m_btOk = new QPushButton("Ok", this);
  m_labActions = new QLabel("Actions:", this);
  m_filesListModel = new QStringListModel(this);
  m_filesListView = new QListView(this);
  m_comboActions = new QComboBox(this);
  m_buttonsWidget = new QWidget(this);
  m_actionsLayout = new QVBoxLayout(m_buttonsWidget);
  m_buttonsLayout = new QHBoxLayout();
  m_mainLayout = new QGridLayout(this);
  
  //  this->openFileDialog->setIncludeFiles(includeFiles);
  //  this->openFileDialog->setExtensions(extensions);
  //  this->openFileDialog->setIncludeNoExtension(includeNoExtension);
  
  m_comboActionItems[ ITEM_CLEAR_LIST ] = "Clear list";
  m_comboActionItems[ ITEM_CLEAR_SELECTION ] = "Clear selection";
  m_comboActionItems[ ITEM_INVERT_SELECTION ] = "Invert selection";
  
  m_filesListView->setModel(m_filesListModel);
  m_filesListView->setSelectionMode(QAbstractItemView::ExtendedSelection);
  
  m_actionsLayout->addWidget(m_labActions);
  
  m_buttonsLayout->addWidget(m_comboActions, 0, Qt::AlignTop);
  m_buttonsLayout->addWidget(m_btOk, 0, Qt::AlignTop);
  
  m_actionsLayout->addWidget(m_labActions, 0);
  m_actionsLayout->addLayout(m_buttonsLayout, 1);
  
  m_comboActions->addItem(m_comboActionItems[ ITEM_ADD_ITEMS ]);
  m_comboActions->addItem(m_comboActionItems[ ITEM_REMOVE_ITEMS ]);
  m_comboActions->addItem(m_comboActionItems[ ITEM_CLEAR_LIST ]);
  m_comboActions->addItem(m_comboActionItems[ ITEM_CLEAR_SELECTION ]);
  m_comboActions->addItem(m_comboActionItems[ ITEM_INVERT_SELECTION ]);
  
  m_mainLayout->addWidget(m_filesListView, 0, 0);
  m_mainLayout->addWidget(m_buttonsWidget, 0, 1);
  
  m_mainLayout->setColumnStretch(0, 1);
  m_mainLayout->setColumnStretch(1, 0);
  
  m_filesListView->setAlternatingRowColors(true);
  
  this->setButtonNames("files");
  
  m_filesListView->setFixedHeight(m_btOk->height() * 4);
  
  connect(m_btOk, SIGNAL(clicked()), this, SLOT(btOkClicked()));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FilesPanel::~FilesPanel()
{
  delete m_filesListView;
  delete m_filesListModel;
  delete m_mainLayout;
  delete m_buttonsLayout;
  delete m_buttonsWidget;
  //delete this->openFileDialog;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QStringList FilesPanel::getFilesList() const
{
  return m_filesListModel->stringList(); 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FilesPanel::setFilesList(const QStringList & filesList)
{
  m_filesListModel->setStringList(filesList);
} 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FilesPanel::addFile()
{
  QStringList files;// = this->openFileDialog->showMultipleSelect();
  QStringList currentFilesList = this->getFilesList();
  
  QStringList::iterator it = files.begin();
  
  while(it != files.end())
  {
    QString filename = *it; 
    
    // if the file is not already in the list
    if(!currentFilesList.contains(filename))
      currentFilesList << filename;
    
    it++;
  }
  
  m_filesListModel->setStringList(currentFilesList);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FilesPanel::removeFile()
{
  QModelIndexList selectedItems = m_filesListView->selectionModel()->
  selectedIndexes();
  
  for(int i = selectedItems.size() - 1 ; i >= 0 ; i--)
  {
    QModelIndex index = selectedItems.at(i);
    
    m_filesListModel->removeRow(index.row(), index.parent());
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FilesPanel::clearList()
{
  m_filesListModel->setStringList(QStringList());
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FilesPanel::setButtonNames(const QString & name)
{
  m_comboActions->setItemText(0, QString("Add %1").arg(name));
  m_comboActions->setItemText(1, QString("Remove %1").arg(name));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void FilesPanel::invertSelection()
{
  QItemSelectionModel * selectionModel = m_filesListView->selectionModel();
  
  for(int i = 0 ; i < m_filesListModel->rowCount() ; i++)
  {
    QModelIndex index = m_filesListModel->index(i, 0);
    
    if(index.isValid())
      selectionModel->select(index, QItemSelectionModel::Toggle);
  }
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void FilesPanel::btOkClicked()
{
  QString selectedText = m_comboActions->currentText();
  int currentIndex = m_comboActionItems.key(selectedText, -1);
  
  if(currentIndex != -1)
  {
    switch(currentIndex)
    {
      case ITEM_CLEAR_LIST:
        this->clearList();
        break;
        
      case ITEM_CLEAR_SELECTION:
        m_filesListView->clearSelection();
        break;
        
      case ITEM_INVERT_SELECTION:
        this->invertSelection();
        break;
    }
  }
  else
  {
    if(selectedText.startsWith("Add"))
      this->addFile();
    
    else if(selectedText.startsWith("Remove"))
      this->removeFile();
  }
}
