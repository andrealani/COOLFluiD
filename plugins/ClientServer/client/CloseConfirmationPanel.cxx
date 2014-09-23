#include <QtGui>

#include "ClientServer/client/CloseConfirmationPanel.hh"

using namespace COOLFluiD::client;

CloseConfirmationPanel::CloseConfirmationPanel(const QString title, 
                                               QDialog * parent)
: QWidget(parent)
{
  m_labTitle = new QLabel(QString("<a href=\".\">%1</a>").arg(title), this);
  m_labText = new QLabel("", this);
  m_buttonsWidget = new QWidget(this);
  m_layout = new QVBoxLayout(this);
  m_buttonsLayout = new QGridLayout(m_buttonsWidget);
  m_buttons = new QDialogButtonBox(this);
  
  m_besideButtonsWidget = NULL;
  
  m_buttons->addButton(QDialogButtonBox::Help);
  
  connect(m_buttons, SIGNAL(helpRequested()), this, SLOT(help()));
  connect(m_labTitle, SIGNAL(linkActivated(const QString &)), 
          this, SLOT(titleClicked(const QString &)));
  
  m_hidden = false;
  
  m_labText->setWordWrap(true);
  m_labTitle->setTextInteractionFlags(Qt::LinksAccessibleByMouse);
  
  m_layout->addWidget(m_labTitle);
  m_layout->addWidget(m_labText);
  
  m_buttonsLayout->addWidget(m_buttons, 0, 1);
  
  m_layout->addWidget(m_buttonsWidget);
  
  this->setLayout(m_layout);
  
  this->setFixedWidth(500); 
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CloseConfirmationPanel::~CloseConfirmationPanel()
{
  /// @todo delete m_buttons
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationPanel::setText(const QString & text)
{
  m_labText->setText(text);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationPanel::setHelp(const QString & helpText)
{
  if(helpText.isEmpty())
    m_helpText = "No help available.";
  
  else
    m_helpText = helpText;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationPanel::addButton(const QString & text, 
                                       QDialogButtonBox::ButtonRole role,
                                       const char * slot)
{
  QPushButton * btn = m_buttons->addButton(text, role);
  connect(btn, SIGNAL(clicked()), this, slot);
  
  m_buttonsVector.push_back(btn);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationPanel::addWidget(QWidget * widget)
{
  if(widget != NULL)
    m_layout->insertWidget(2, widget);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationPanel::hideComponents(bool hide)
{
  m_hidden = hide;
  
  m_buttons->setHidden(hide);
  m_labText->setHidden(hide);
  
  this->adjustSize();
  emit resized();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationPanel::setWidgetBesideButtons(QWidget * widget)
{
  if(m_besideButtonsWidget == NULL && widget != NULL)
    m_buttonsLayout->addWidget(widget, 0, 0);
}

/****************************************************************************
 
 SLOTS
 
 ****************************************************************************/

void CloseConfirmationPanel::help()
{
  emit showHelp(m_helpText);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CloseConfirmationPanel::titleClicked(const QString & link)
{
  this->hideComponents(!m_hidden);
}
