#include <QtGui>


#include "ClientServer/client/CloseConfirmationInfos.hh"
#include "ClientServer/client/RemoteSaveFile.hh"
#include "ClientServer/client/SelectFileDialog.hh"

#include "ClientServer/client/SaveConfirmationPanel.hh"

using namespace COOLFluiD::client;

SaveConfirmationPanel::SaveConfirmationPanel(QDialog * parent, bool becauseCommit)
: CloseConfirmationPanel("Save modified configuration", parent)
{
  m_comboBox = new QComboBox(this);
  m_labStatus = new QLabel("", this);
  m_checkSaveLocallyOnError = new QCheckBox("Save locally on error", this);
  m_fileNameWidget = new QWidget(this);
  m_widgetLayout = new QHBoxLayout(m_fileNameWidget);
  m_editFileName = new QLineEdit(m_fileNameWidget);
  
  if(becauseCommit)
  {
    this->setText("Some options have been modified but not committed. Committing "
                  "them implies that the current configuration will be modified. "
                  "<i>What do you want to do?</i>");
  }
  else
  {
    this->setText("The current configuration has been modified but not saved. "
                  "<i>What do you want to do?</i>");
  }
  
  QString help = "You can choose to save the configuration locally (this "
  "computer) or remotely (the computer on which the server is running). If "
  "you choose to save remotely, you have the ability to check \"Save locally "
  "on error\" checkbox."
  "\nChecking it means that the file has to be saved locally if it "
  "couldn't be saved remotely. In that case, the file will be saved in \"%1\" "
  "with the same name; a number will be added to the filename if needed to "
  "avoid overwriting an existing file.\n";
  
  if(becauseCommit)
  {
    help.prepend("The current configuration has not been modified, but it will "
                 "be if you choose to commit modified options. Thus, you have here the "
                 "possibility to save the configuration, if needed.\n");
  }
  
  m_rollingBack = false;
  m_previousIndex = 0;
  
  m_labStatus->setBuddy(m_editFileName);
  m_editFileName->setReadOnly(true);
  m_editFileName->setVisible(false);
  
  this->setHelp(help.arg(QDir::home().path()));
  
  m_comboBox->addItem("--> [ Select an item ] <--");
  m_comboBox->addItem("Don't save");
  m_comboBox->addItem("Save locally");
  m_comboBox->addItem("Save remotely");
  
  m_checkSaveLocallyOnError->setEnabled(false);
  
  m_widgetLayout->addWidget(m_labStatus);
  m_widgetLayout->addWidget(m_editFileName);
  
  this->addWidget(m_fileNameWidget);
  this->addWidget(m_comboBox);
  
  this->setWidgetBesideButtons(m_checkSaveLocallyOnError);
  
  connect(m_comboBox, SIGNAL(currentIndexChanged(int)), 
          this, SLOT(currentIndexChanged(int)));
  
  this->hideComponents(true);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SaveConfirmationPanel::~SaveConfirmationPanel()
{
  delete m_comboBox;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool SaveConfirmationPanel::isAccepted() const
{
  return m_comboBox->currentIndex() != 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SaveConfirmationPanel::getData(CloseConfirmationInfos & infos) const
{
  int index = m_comboBox->currentIndex();
  
  if(this->isAccepted())
  {
    infos.saveLocally = index == 2;
    infos.saveLocallyOnError = m_checkSaveLocallyOnError->isEnabled() && 
    m_checkSaveLocallyOnError->isChecked();
    
    infos.filename = index <= 1 ? "" : m_filename; 
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SaveConfirmationPanel::setData(const CloseConfirmationInfos & infos)
{
  /// @todo write the code
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SaveConfirmationPanel::hideComponents(bool hide)
{
  m_comboBox->setHidden(hide);
  m_checkSaveLocallyOnError->setHidden(hide);
  m_fileNameWidget->setHidden(hide);
  
  CloseConfirmationPanel::hideComponents(hide);
}

/****************************************************************************
 
 PRIVATE METHODS
 
 ****************************************************************************/

QString SaveConfirmationPanel::saveLocally()
{
  SelectFileDialog sfd;
  sfd.addFileType("CFcase", "CFcase");
  sfd.addFileType("XML", "xml");
  return sfd.show(QFileDialog::AcceptSave);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString SaveConfirmationPanel::saveRemotely()
{
  //  RemoteSaveFile rsf;
  //  rsf.setExtensions(QStringList() << "CFcase" << "xml");
  //  return rsf.show();
  
  return "";
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void SaveConfirmationPanel::currentIndexChanged(int index)
{
  m_checkSaveLocallyOnError->setEnabled(false);
  
  if(index == 0)
  {
    m_filename.clear();
    m_labStatus->setText("");
    m_previousIndex = 0;
    m_editFileName->setVisible(false);
  }
  
  else if(index == 1)
  {
    m_filename.clear();
    m_labStatus->setText("Configuration will not be saved!");
    m_previousIndex = 1;
    m_editFileName->setVisible(false);
  }
  
  else if(index == 2)
  {
    QString filename = this->saveLocally();
    // repaint the dialog that contains this panel (this was not done automatically)
    this->parentWidget()->repaint();
    
    if(filename.isEmpty())
    {
      m_rollingBack = true;
      m_comboBox->setCurrentIndex(m_previousIndex);
    }
    else
    {
      m_filename = filename;
      m_labStatus->setText("Save to (locally):");
      m_editFileName->setText(filename);
      m_editFileName->setVisible(true);
      m_previousIndex = 2;
    }
  }
  
  else if(index == 3)
  {
    QString filename = this->saveRemotely();
    // repaint the dialog that contains this panel (this was not done automatically)
    this->parentWidget()->repaint(); 
    
    if(filename.isEmpty())
    {
      m_rollingBack = true;
      m_comboBox->setCurrentIndex(m_previousIndex);
    }
    else
    {
      m_filename = filename;
      m_labStatus->setText("Save to (remotely):");
      m_editFileName->setText(filename);
      m_editFileName->setVisible(true);
      m_previousIndex = 3;
    }
    
    m_checkSaveLocallyOnError->setEnabled(true);
  }
  
  this->adjustSize();
  emit resized();
}
