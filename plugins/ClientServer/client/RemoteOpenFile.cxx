#include <QtGui>

#include "ClientServer/client/FilesListItem.hh"
#include "ClientServer/client/RemoteOpenFile.hh"

using namespace COOLFluiD::client;

RemoteOpenFile::RemoteOpenFile(const QModelIndex & index, QMainWindow * parent)
: RemoteFSBrowser(index, parent)
{
  this->setIncludeFiles(true);
  this->setIncludeNoExtension(true);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

RemoteOpenFile::~RemoteOpenFile()
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ValidationPolicy RemoteOpenFile::isAcceptable(const QString & name, bool isDir) 
{
  if(isDir)
    return POLICY_ENTER_DIRECTORY;
  
  m_fileList << name;
  return POLICY_VALID;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ValidationPolicy RemoteOpenFile::isAcceptable(const QStringList & names) 
{
  QStringList::const_iterator it = names.begin();
  ValidationPolicy validation = POLICY_VALID;
  
  while(it != names.end() && validation == POLICY_VALID)
  {
    QString item = *it;
    
    if(this->isDirectory(item) && names.size() > 1)
    {
      this->showError("Directories are not allowed in multiple selection.");
      validation = POLICY_NOT_VALID;
      m_fileList.clear();
    }
    
    else if(this->isDirectory(item) && names.size() == 1)
      validation = POLICY_ENTER_DIRECTORY;
    
    else
      m_fileList << item;
    
    it++;
  }
  
  
  return validation;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString RemoteOpenFile::getSelectedFile() const
{
  if(!m_fileList.isEmpty())
    return m_fileList.at(0);
  
  return QString();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QStringList RemoteOpenFile::getSelectedFileList() const
{
  return m_fileList;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RemoteOpenFile::reinitValues()
{
  m_fileList.clear();
}
