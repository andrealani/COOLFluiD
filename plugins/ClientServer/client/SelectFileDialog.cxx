#include <QtGui>

#include "ClientServer/client/SelectFileDialog.hh"

using namespace COOLFluiD::client;

SelectFileDialog::SelectFileDialog(QWidget * parent)
: QFileDialog(parent)
{
  this->setDirectory(QDir::home());
  this->setViewMode(QFileDialog::Detail);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SelectFileDialog::~SelectFileDialog()
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString SelectFileDialog::show(AcceptMode mode)
{
  QStringList filters;
  QHash<QString, QString>::iterator it = m_extensions.begin();
  
  while(it != m_extensions.end())
  {
    QString value = it.value().replace(" ", " *."); // set "*." before each extension
    QString key = it.key();
    
    if(!value.isEmpty())
      filters << QString("%1 files (*.%2)").arg(key).arg(value);
    else
      filters << QString("%1 files (*)").arg(key);
    
    it++; 
  }
  
  this->setNameFilters(filters);
  this->setAcceptMode(mode); // dialog title is automatically changed
  
  // if "open", the file must exist
  if(mode == QFileDialog::AcceptOpen)
    this->setFileMode(QFileDialog::ExistingFile);
  else 
    this->setFileMode(QFileDialog::AnyFile);
  
  if(this->exec() == QDialog::Accepted)
    return m_filename;
  else
    return QString();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SelectFileDialog::addFileType(const QString & type, 
                                   const QString & extensions)
{
  QString exts = this->unique(extensions);
  QString fileType = type.trimmed();
  
  if(!fileType.isEmpty() && !exts.isEmpty())
  {
    exts.remove(QRegExp("[\\*\\?]"));
    
    m_extensions[fileType] = exts;
  }
  
  else if(fileType.isEmpty() && exts.isEmpty())
    m_extensions["All"] = "";
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SelectFileDialog::removeFileType(const QString & type)
{
  m_extensions.remove(type);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SelectFileDialog::appendExtensionToFileType(const QString & type, 
                                                 const QString & extensions)
{
  QString trimmedType = type.trimmed().remove(QRegExp("[\\*\\?]"));
  if(!trimmedType.isEmpty())
  {
    QString newExts = this->unique(m_extensions[type] + " " + extensions);
    m_extensions[trimmedType] = newExts;
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SelectFileDialog::clearFileTypes()
{
  m_extensions.clear();
}

/****************************************************************************
 
 SLOT
 
 ****************************************************************************/

void SelectFileDialog::accept()
{
  QString regexPattern;
  QRegExp regex;
  QStringList files = this->selectedFiles();
  QHash<QString, QString>::iterator it = m_extensions.begin();
  
  // if no file is selected, nothing to do
  if(files.isEmpty())
    return;
  
  // build a regular expression to check if the selected file has one of the
  // listed extensions. i.e. : if listed extension are xml, CFcase and txt, the 
  // resulting regular expression will be : "^.+\\.((xml)|(CFcase)|(txt))$"
  while(it != m_extensions.end())
  {
    if(!regexPattern.isEmpty())
      regexPattern.append(")|(");
    
    regexPattern.append(it.value());
    it++;
  }
  
  if(!regexPattern.isEmpty())
    regexPattern.prepend("^.+\\.((").append("))$");
  
  // if no extension is listed, all extensions are accepted...call base method
  else
    QFileDialog::accept();
  
  regex.setPattern(regexPattern);
  // multiple file selection not allowed, if a file is selected, it is always
  // the first string of the string list
  m_filename = files.at(0); 
  
  // if this is a "Save dialog" and the selection is confirmed, call base method
  if(this->acceptMode() == QFileDialog::AcceptSave && this->fixExtensionSave(regex))
    QFileDialog::accept(); // asks to confirm if overwrite an existing file
  
  else if(this->acceptMode() == QFileDialog::AcceptOpen)
  {
    this->fixExtensionOpen();
    QFileDialog::accept(); // fails if the selected file does not exist
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool SelectFileDialog::fixExtensionSave(const QRegExp & extensions)
{
  bool accept = true;
  int extPosition = m_filename.length() - m_filename.lastIndexOf('.') - 1;
  QString fileExt = m_filename.right(extPosition);
  QString selectedType = this->selectedNameFilter();
  QRegExp regexType = QRegExp(" files (*)", Qt::CaseSensitive, QRegExp::Wildcard);
  QString key = selectedType.remove(regexType);
  QString keyForExtension = m_extensions.key(fileExt);
  
  // if the file has one of the extensions listed in the regular expression
  if(m_filename.contains(extensions)) 
  {
    // if the extension if the one of the selected type (i.e. : "xml" extension 
    // but type is "CFcase"), ask the user to choose to keep the selected type,
    // switch the type of the extension or cancel
    if(m_extensions[key] != fileExt && !keyForExtension.isEmpty())
    {
      QMessageBox confirmDialog(this);
      QString text = "You selected %1 file type but the name you entered has %2 "
      "extension.";
      QString tmpRegex = QString("\\.((%1)|(%2))$").arg(m_extensions[key]).arg(fileExt);
      QPushButton * btKeep;
      QPushButton * btSwitch;
      
      // configure the dialog
      confirmDialog.setWindowTitle("File type conflict");
      
      confirmDialog.setText(text.arg(key).arg(fileExt));
      
      confirmDialog.setInformativeText("What do you want to do?");
      
      btKeep = confirmDialog.addButton(QString("Keep %1 type").arg(key), 
                                       QMessageBox::YesRole);
      btSwitch = confirmDialog.addButton(QString("Switch to %1 type").arg(keyForExtension), 
                                         QMessageBox::NoRole);
      
      confirmDialog.addButton(QMessageBox::Cancel);
      
      // if user cancels his entry
      if(confirmDialog.exec() != QMessageBox::Cancel)
      {
        /// @todo don't remove extension if btSwitch
        
        // remove the extension
        m_filename.remove(QRegExp(tmpRegex));
        
        // if user wants to keep the type, set the correct extension
        if(confirmDialog.clickedButton() == btKeep)
          m_filename.append("." + m_extensions[key]);
        
        else if(confirmDialog.clickedButton() == btSwitch)
        {
          m_filename.append("." + fileExt);
          // select the correct type
          this->selectNameFilter(QString("%1 files (*.%2)").arg(keyForExtension)
                                 .arg(m_extensions[keyForExtension]));
        }
      } // END "if(confirmDialog.exec() != QMessageBox::Cancel)"
      else
        accept = false;
      
    } // END "if(this->extensions[key] != fileExt && !keyForExtension.isEmpty())"
    else if(keyForExtension.isEmpty())
      m_filename.append("." + m_extensions[key]);
    
    QMessageBox::information(this, NULL, m_filename);
    
  } // END "if(this->filename.contains(extensions))"
  else
    m_filename.append("." + m_extensions[key]);
  
  return accept;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SelectFileDialog::fixExtensionOpen()
{
  QString selectedType = this->selectedNameFilter();
  QRegExp regexType = QRegExp(" files (*)", Qt::CaseSensitive, QRegExp::Wildcard);
  QString key = selectedType.remove(regexType);
  
  // if the file does not have the right extension, remove the wrong one (if any)
  // and append the correct one
  if(!m_filename.endsWith(m_extensions[key])) 
  {
    int extPosition = m_filename.length() - m_filename.lastIndexOf('.') - 1;
    QString fileExt = m_filename.right(extPosition);
    
    if(m_extensions[key] != fileExt)
      m_filename.append("." + m_extensions[key]);
  } 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString SelectFileDialog::unique(const QString & string)
{
  QStringList list = string.trimmed().split(" ", QString::SkipEmptyParts);
  list.removeDuplicates();
  return list.join(" ");
}
