#include <QtGui>

#include "ClientServer/client/CloseConfirmationInfos.hh"
#include "ClientServer/client/CommitConfirmationPanel.hh"

using namespace COOLFluiD::client;

CommitConfirmationPanel::CommitConfirmationPanel(QDialog * parent)
: CloseConfirmationPanel("Commit modified options", parent)
{
  this->comboBox = new QComboBox(this);
  
  this->addButton("Details", QDialogButtonBox::ActionRole, SLOT(showDetails()));
  
  this->setText("Options were modified but not committed. <i>What do you want "
                "to do?</i>");
  
  this->setHelp("You can click on \"Details\" button to that see options that "
                "have been modified (with new and old value) and the ones that have been "
                "added.");
  
  this->comboBox->addItem("--> [ Select an item ] <--");
  this->comboBox->addItem("Don't commit");
  this->comboBox->addItem("Commit");
  
  this->addWidget(this->comboBox);
  
  this->hideComponents(true);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CommitConfirmationPanel::~CommitConfirmationPanel()
{
  delete this->comboBox;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitConfirmationPanel::isAccepted() const
{
  return this->comboBox->currentIndex() != 0;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitConfirmationPanel::getData(CloseConfirmationInfos & infos) const
{
  
  if(this->isAccepted())
    infos.commitRequested = this->comboBox->currentIndex() == 2;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitConfirmationPanel::setData(const CloseConfirmationInfos & infos)
{
  QString separator = QString("\n").rightJustified(30, '+');
  
  this->details.clear();
  
  if(infos.commitRequested)
    this->comboBox->setCurrentIndex(2);
  
  else
    this->comboBox->setCurrentIndex(0);
  
  this->details = infos.commitDetails.toString(); 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitConfirmationPanel::hideComponents(bool hide)
{
  this->comboBox->setHidden(hide);
  CloseConfirmationPanel::hideComponents(hide);
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void CommitConfirmationPanel::showDetails()
{
  emit showHelp(this->details);
}