#include <QtGui>

#include "ClientServer/client/ClientKernel.hh"
#include "ClientServer/client/TypeAndNameDialog.hh"
#include "ClientServer/client/RemoteSaveFile.hh"

using namespace COOLFluiD::client;

RemoteSaveFile::RemoteSaveFile(const QModelIndex & index, QMainWindow * parent)
: RemoteFSBrowser(index, parent)
{
  this->setIncludeFiles(true);
  this->setIncludeNoExtension(false);
  
  this->setWindowTitle("Save configuration");
  
  this->allowMultipleSelect = false;
  
  m_btFileName = this->addButton("Set file name", 
                                       QDialogButtonBox::ActionRole);
  
  m_btNewDirectory = this->addButton("New directory", 
                                           QDialogButtonBox::ActionRole);
  
  m_fileNameDialog = new TypeAndNameDialog("File name", "File extension", 
                                                 (QWidget *) this);
  
  connect(m_btFileName, SIGNAL(clicked()), 
          this, SLOT(btFileNameClick()));
  
  connect(m_btNewDirectory, SIGNAL(clicked()), 
          this, SLOT(btNewDirectoryClicked()));
  
  this->allowModifyBools = false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

RemoteSaveFile::~RemoteSaveFile()
{
  delete m_fileNameDialog;
  delete m_btFileName;
  delete m_btNewDirectory;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RemoteSaveFile::btFileNameClick()
{
  QString selectedExtension;
  QString name;
  bool overwrite;
  
  do
  {
    overwrite = true;
    name = m_fileNameDialog->show(this->getExtensions(), selectedExtension);
    
    if(!name.isEmpty())
    {
      if(!name.endsWith(selectedExtension))
        name.append(".").append(selectedExtension);
      
      this->setStatus(QString("Current file name: \"%1\"").arg(name));
      
      if(this->itemExists(name))
      {
        int answer;
        QString path = this->getCurrentPath();
        QString message = "The file '%1' already exists. Are you sure to "
        "overwrite this file?";
        
        this->assemblePath(path, name);
        
        answer = QMessageBox::question(this, "Confirmation", message.arg(path),
                                       QMessageBox::Yes | QMessageBox::No);
        
        overwrite = answer == QMessageBox::Yes;
      }
    }
  } while(!overwrite); // we continue while user says to not overwrite
  
  m_fileName = name;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RemoteSaveFile::btNewDirectoryClicked()
{
  QString dirName;
  bool ok;
  
  dirName = QInputDialog::getText(this, "New directory", "Enter the name of the "
                                  "new directory", QLineEdit::Normal, "", &ok);
  
  if(ok)
  {
    if(!dirName.isEmpty())
      this->clientKernel->createDir(this->index, this->getCurrentPath(), dirName);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ValidationPolicy RemoteSaveFile::isAcceptable(const QString & name, bool isDir)
{
  ValidationPolicy validation = POLICY_NOT_VALID;
  
  if(isDir && m_fileName.isEmpty())
  {
    this->showError(QString("You must select a file or enter a new file name "
                            "using '%1'").arg(m_btFileName->text()));
  }
  
  else if(isDir)
  {
    // this->selectedFile is empty, because reinitValues() was called
    m_selectedFile.append(name);
    this->assemblePath(m_selectedFile, m_fileName);
    validation = POLICY_VALID;
  }
  
  else if(!isDir)
  {
    int answer = QMessageBox::question(this, "Confirmation", QString("The file "
                                                                     "'%1' already exists. Are you sure you want to overwrite this "
                                                                     "file?").arg(name), QMessageBox::Yes | QMessageBox::No);
    
    if(answer == QMessageBox::Yes)
    {
      validation = POLICY_VALID;
      m_selectedFile = name;
    }
  }
  
  return validation;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RemoteSaveFile::reinitValues()
{
  m_selectedFile.clear();
  m_fileName.clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString RemoteSaveFile::getSelectedFile() const
{
  return m_selectedFile;
}
