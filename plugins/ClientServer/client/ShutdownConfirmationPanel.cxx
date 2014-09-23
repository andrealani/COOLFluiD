#include <QtGui>

#include "ClientServer/client/CloseConfirmationInfos.hh"
#include "ClientServer/client/ShutdownConfirmationPanel.hh"

using namespace COOLFluiD::client;

ShutdownConfirmationPanel::ShutdownConfirmationPanel(QDialog * parent)
: CloseConfirmationPanel("Shutdown the server", parent)
{
  m_comboBox = new QComboBox(this);
  
  this->setText("The server is still running. <i>What do you want to do?</i>");
  
  this->setHelp("If you choose to shutdown the server while a simulation is "
                "running, the simulation will be stopped and data may be lost.");
  
  m_comboBox->addItem("--> [ Select an item ] <--");
  m_comboBox->addItem("Keep the server running");
  m_comboBox->addItem("Shutdown the server");
  
  this->addWidget(m_comboBox);
  
  this->hideComponents(true);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ShutdownConfirmationPanel::~ShutdownConfirmationPanel()
{
  delete m_comboBox;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ShutdownConfirmationPanel::isAccepted() const
{
  return m_comboBox->currentIndex() != 0;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ShutdownConfirmationPanel::getData(CloseConfirmationInfos & infos) const
{
  if(this->isAccepted())
    infos.shutdownServerRequested = m_comboBox->currentIndex() == 2;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ShutdownConfirmationPanel::setData(const CloseConfirmationInfos & infos)
{
  if(infos.shutdownServerRequested)
    m_comboBox->setCurrentIndex(2);
  
  else
    m_comboBox->setCurrentIndex(0);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ShutdownConfirmationPanel::hideComponents(bool hide)
{
  m_comboBox->setHidden(hide);
  CloseConfirmationPanel::hideComponents(hide);
}
