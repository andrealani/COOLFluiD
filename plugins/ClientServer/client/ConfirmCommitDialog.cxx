#include <QtGui>

#include "ClientServer/client/CommitDetails.hh"

#include "ClientServer/client/ConfirmCommitDialog.hh"

using namespace COOLFluiD::client;

ConfirmCommitDialog::ConfirmCommitDialog(QWidget * parent)
: QDialog(parent)
{
  this->setWindowTitle("Commit confirm");
  
  m_labText = new QLabel("Options have been modified but were not comitted.<br>"
                               "Click on \"<i>Details</i>\" to see what "
                               "modifications have been done.", this);
  
  
  m_mainLayout = new QVBoxLayout(this);
  
  m_buttonBox = new QDialogButtonBox(this);
  
  m_detailsView = new QTableView(this);
  
  this->createButton("Cancel", CANCEL, QDialogButtonBox::RejectRole);
  this->createButton("Commit", COMMIT, QDialogButtonBox::YesRole);
  this->createButton("Do not commit", DONT_COMMIT, QDialogButtonBox::NoRole);
  
  m_mainLayout->addWidget(m_labText);
  m_mainLayout->addWidget(m_detailsView);
  m_mainLayout->addWidget(m_buttonBox);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ConfirmCommitDialog::~ConfirmCommitDialog()
{
  
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ConfirmCommitDialog::CommitConfirmation ConfirmCommitDialog::show(CommitDetails & commitDetails)
{
  m_answer = CANCEL;
  
  if(commitDetails.isEmpty())
  { 
    m_detailsView->setModel(&commitDetails);
    this->exec();
  }
  
  return m_answer;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ConfirmCommitDialog::buttonClicked()
{
  QPushButton * button = static_cast<QPushButton *> (sender());
  
  if(button != NULL)
    m_answer = m_buttons.key(button);
  else 
    m_answer = CANCEL;
  
  this->hide();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ConfirmCommitDialog::createButton(const QString & text, 
                                       CommitConfirmation commConf, 
                                       QDialogButtonBox::ButtonRole role)
{
  QPushButton * button = m_buttonBox->addButton(text, role);
  connect(button, SIGNAL(clicked()), this, SLOT(buttonClicked()));
  m_buttons[commConf] = button;
}
