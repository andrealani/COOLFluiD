#include <QtCore>
#include <QtGui>

#include "ClientServer/client/RemoteFSBrowser.hh"
#include "ClientServer/client/LibrariesPanel.hh"

using namespace COOLFluiD::client;

LibrariesPanel::LibrariesPanel(QWidget * parent) 
: FilesPanel(true, QStringList() << "so", false, parent)
{ 
  this->setButtonNames("libraries");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

LibrariesPanel::~LibrariesPanel()
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QStringList LibrariesPanel::getFilesList() const
{
  QStringList libraries;
  QStringList files = FilesPanel::getFilesList();
  QStringList::iterator it = files.begin();
  
  while(it != files.end())
  {
    QFileInfo path(*it);
    /// @todo make platform independant
    // get the file name and remove extension
    // note : QFileInfo::baseName() returns directly the file name without
    // extension but I noticed that it may remove a part of the file name too 
    // (if it contains the version of the library). 
    // e.g. : for "/lib/libpthread-2.8.so", baseName() returns "libpthread-2"
    libraries << path.fileName().remove(QRegExp("\\.so$"));
    it++;
  }
  
  return libraries; 
}
