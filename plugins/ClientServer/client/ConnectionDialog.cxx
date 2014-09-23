#include <QtGui>

#include "ClientServer/client/ConnectionDialog.hh"
#include "ClientServer/treeview/TSshInformation.hh"

using namespace COOLFluiD::client;
using namespace COOLFluiD::treeview;

ConnectionDialog::ConnectionDialog(QMainWindow * parent)
: QDialog(parent)
{
  QString username;
  QRegExp regex("^USER=");
  QStringList environment = QProcess::systemEnvironment().filter(regex);
  
  if(environment.size() == 1)
    username = environment.at(0);
  
  this->setWindowTitle("Connect to server");
  
  // create the components
  m_labHostname = new QLabel("Hostname:");
  m_labUsername = new QLabel("Username:");
  m_labPortNumber = new QLabel("Port number:");
  
  m_editHostname = new QLineEdit(this);
  m_editUsername = new QLineEdit(this);
  m_spinPortNumber = new QSpinBox(this);
  
  m_infosLayout = new QHBoxLayout();
  
  m_chkLaunchServer = new QCheckBox("Start a new server instance", this);
  
  m_layout = new QFormLayout(this);
  m_buttons = new QDialogButtonBox(QDialogButtonBox::Ok
                                         | QDialogButtonBox::Cancel);
  
  // the dialog is modal
  this->setModal(true);
  
  m_spinPortNumber->setMinimum(49150);
  m_spinPortNumber->setMaximum(65535);
  
  m_editHostname->setText("localhost");
  m_editUsername->setText(username.remove("USER="));
  m_spinPortNumber->setValue(62784);
  
  // add the components to the m_layout
  m_infosLayout->addWidget(m_labHostname);
  m_infosLayout->addWidget(m_editHostname);
  m_infosLayout->addWidget(m_labPortNumber);
  m_infosLayout->addWidget(m_spinPortNumber); 
  
  this->chkLaunchServerChecked(m_chkLaunchServer->checkState());
  
  m_layout->addRow(m_infosLayout);
  m_layout->addRow(m_chkLaunchServer);
  m_layout->addRow(m_labUsername, m_editUsername);
  m_layout->addRow(m_buttons);
  
  // add the m_layout to the dialog
  this->setLayout(m_layout);
  
  // connect useful signals to slots
  connect(m_buttons, SIGNAL(accepted()), this, SLOT(btOkClicked()));
  connect(m_buttons, SIGNAL(rejected()), this, SLOT(btCancelClicked()));
  connect(m_chkLaunchServer, SIGNAL(stateChanged(int)),
          this, SLOT(chkLaunchServerChecked(int)));
  
  //  m_minSize = this->size();
  m_layout->setSizeConstraint(QLayout::SetMaximumSize);
  //  this->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ConnectionDialog::~ConnectionDialog()
{
  delete m_buttons;
  delete m_chkLaunchServer;
  delete m_editUsername;
  delete m_editHostname;
  delete m_infosLayout;
  delete m_labHostname;
  delete m_labPortNumber;
  delete m_labUsername;
  delete m_layout;
  delete m_spinPortNumber;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ConnectionDialog::show(bool hidePort, TSshInformation & sshInfos)
{
  m_okClicked = false;
  
  m_labPortNumber->setVisible(!hidePort);
  m_spinPortNumber->setVisible(!hidePort);
  this->adjustSize();
  
  this->exec();
  
  if(m_okClicked)
  {
    sshInfos.m_hostname = m_editHostname->text();
    sshInfos.username = m_editUsername->text();
    sshInfos.launchServer = m_chkLaunchServer->isChecked();
    sshInfos.port = m_spinPortNumber->value();
  }
  
  return m_okClicked;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ConnectionDialog::setSshInfos(const TSshInformation & sshInfos)
{
  Qt::CheckState checked = sshInfos.launchServer ? Qt::Checked : Qt::Unchecked;
  
  m_editHostname->setText(sshInfos.m_hostname);
  m_editUsername->setText(sshInfos.username);
  m_chkLaunchServer->setCheckState(checked);
  m_spinPortNumber->setValue(sshInfos.port);
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void ConnectionDialog::btOkClicked()
{
  m_okClicked = true;
  this->setVisible(false);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ConnectionDialog::btCancelClicked()
{
  this->setVisible(false);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ConnectionDialog::chkLaunchServerChecked(int state)
{
  m_editUsername->setVisible(state == Qt::Checked);
  m_labUsername->setVisible(state == Qt::Checked);
  this->adjustSize();
}
