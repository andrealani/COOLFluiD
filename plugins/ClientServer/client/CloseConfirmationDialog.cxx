#include <QtGui>

#include "ClientServer/client/CommitConfirmationPanel.hh"
#include "ClientServer/client/CloseConfirmationInfos.hh"
#include "ClientServer/client/CloseConfirmationPanel.hh"
#include "ClientServer/client/SaveConfirmationPanel.hh"
#include "ClientServer/client/ShutdownConfirmationPanel.hh"

#include "ClientServer/client/CloseConfirmationDialog.hh"

using namespace COOLFluiD::client;

CloseConfirmationDialog::CloseConfirmationDialog(QMainWindow * parent)
: QDialog(parent)
{
  QString instr = "At least one element needs your attention before you close "
  "the application. Please, answer to <b>all</b> questions below before "
  "continue. Click on a question to show/hide details about it. You can "
  "click on the appropriate \"Help\" button to see some help about an "
  "element.";
  
  this->setWindowTitle("Close confirmation");
  
  m_labInstructions = new QLabel(instr, this);
  m_layout = new QVBoxLayout(this);
  m_buttons = new QDialogButtonBox(QDialogButtonBox::Ok
                                         | QDialogButtonBox::Cancel); 
  m_helpDialog = new QDialog(this);
  m_editHelp = new QTextEdit(m_helpDialog);
  
  m_panels[ CLOSE_SAVE_FILE ] = NULL;
  m_panels[ CLOSE_COMMIT ] = NULL;
  m_panels[ CLOSE_SHUT_DOWN ] = NULL;
  
  m_labInstructions->setWordWrap(true);
  m_labInstructions->setTextFormat(Qt::RichText);
  
  m_editHelp->setWordWrapMode(QTextOption::WordWrap);
  m_editHelp->setReadOnly(true);
  m_helpDialog->setWindowTitle("Help");
  
  m_okClicked = false;
  m_helpNeverShown = true;
  
  m_layout->addWidget(m_labInstructions, 0, 0); 
  
  connect(m_buttons, SIGNAL(accepted()), this, SLOT(btOkClicked()));
  connect(m_buttons, SIGNAL(rejected()), this, SLOT(btCancelClicked()));
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CloseConfirmationDialog::~CloseConfirmationDialog()
{
  delete m_buttons;
  delete m_labInstructions;
  
  delete m_panels[ CLOSE_SAVE_FILE ];
  delete m_panels[ CLOSE_COMMIT ];
  delete m_panels[ CLOSE_SHUT_DOWN ];
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationDialog::addConfirmation(CloseConfirmationType type,
                                              bool becauseCommit)
{
  if(m_panels[type] == NULL)
  {
    if(type == CLOSE_SHUT_DOWN)
      m_panels[type] = new ShutdownConfirmationPanel(this);
    
    if(type == CLOSE_SAVE_FILE)
      m_panels[type] = new SaveConfirmationPanel(this, becauseCommit);
    
    if(type == CLOSE_COMMIT)
      m_panels[type] = new CommitConfirmationPanel(this);
    
    connect(m_panels[type], SIGNAL(showHelp(const QString &)), 
            this, SLOT(showHelp(const QString &)));
    
    connect(m_panels[type], SIGNAL(resized()), this, SLOT(panelResized()));
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CloseConfirmationDialog::show(CloseConfirmationInfos & infos)
{
  if(m_panels[ CLOSE_COMMIT ] == NULL && 
     m_panels[ CLOSE_SAVE_FILE ] == NULL &&
     m_panels[ CLOSE_SHUT_DOWN ] == NULL)
    return true;
  
  if(m_panels[ CLOSE_COMMIT ] != NULL)
  {
    m_layout->addWidget(m_panels[ CLOSE_COMMIT ]);
    m_panels[ CLOSE_COMMIT ]->setData(infos);
  }
  
  if(m_panels[ CLOSE_SAVE_FILE ] != NULL)
  {
    m_layout->addWidget(m_panels[ CLOSE_SAVE_FILE ]);
    m_panels[ CLOSE_SAVE_FILE ]->setData(infos);
  }
  
  if(m_panels[ CLOSE_SHUT_DOWN ] != NULL)
  {
    m_layout->addWidget(m_panels[ CLOSE_SHUT_DOWN ]);
    m_panels[ CLOSE_SHUT_DOWN ]->setData(infos);
  }
  
  
  m_layout->addWidget(m_buttons);
  
  this->adjustSize();
  
  m_editHelp->setFixedWidth(this->width());
  m_editHelp->setFixedHeight(100);
  m_helpDialog->adjustSize();
  m_helpDialog->setFixedSize(m_helpDialog->size());
  
  this->exec();
  
  if(m_okClicked)
  {
    if(m_panels[ CLOSE_SHUT_DOWN ] != NULL)
      m_panels[ CLOSE_SHUT_DOWN ]->getData(infos);
    
    if(m_panels[ CLOSE_SAVE_FILE ] != NULL)
      m_panels[ CLOSE_SAVE_FILE ]->getData(infos);
    
    if(m_panels[ CLOSE_COMMIT ] != NULL)
      m_panels[ CLOSE_COMMIT ]->getData(infos);
  }
  
  return m_okClicked;
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void CloseConfirmationDialog::btOkClicked()
{
  bool allAccepted = true;
  
  if(m_panels[ CLOSE_COMMIT ] != NULL)
    allAccepted &= m_panels[ CLOSE_COMMIT ]->isAccepted();
  
  if(m_panels[ CLOSE_SAVE_FILE ] != NULL)
    allAccepted &= m_panels[ CLOSE_SAVE_FILE ]->isAccepted();
  
  if(m_panels[ CLOSE_SHUT_DOWN ] != NULL)
    allAccepted &= m_panels[ CLOSE_SHUT_DOWN ]->isAccepted();
  
  if(!allAccepted)
    QMessageBox::critical(this, "Error", "You did not answer all questions.");
  
  else
  {
    m_okClicked = true;
    this->setVisible(false);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationDialog::btCancelClicked()
{
  m_okClicked = false;
  this->setVisible(false);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationDialog::showHelp(const QString & helpText)
{
  if(m_helpNeverShown)
  {
    QPoint position = this->frameGeometry().topLeft(); 
    position.setY(position.y() + this->frameGeometry().height());
    m_helpDialog->move(position);
    m_helpNeverShown = false;
  }
  
  m_editHelp->setText(helpText);
  m_helpDialog->show();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationDialog::panelResized()
{
  this->adjustSize();
}
