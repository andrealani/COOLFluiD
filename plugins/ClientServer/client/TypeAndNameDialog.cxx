#include <QtGui>

#include "ClientServer/client/TypeAndNameDialog.hh"

using namespace COOLFluiD::client;

TypeAndNameDialog::TypeAndNameDialog(const QString & fieldLabelText,
                                     const QString & dropListText, 
                                     QWidget * parent)
: QDialog(parent)
{
  this->setWindowTitle("Add a new child node");
  
  // create the components
  m_labName = new QLabel(fieldLabelText + ":");
  m_labConcreteType = new QLabel(dropListText + ":");
  m_editName = new QLineEdit(this);
  m_cbTypes = new QComboBox(this);
  m_layout = new QFormLayout(this); 
  m_buttons = new QDialogButtonBox(QDialogButtonBox::Ok | 
                                         QDialogButtonBox::Cancel);
  
  // add the components to the m_layout
  m_layout->addRow(m_labName, m_editName);
  m_layout->addRow(m_labConcreteType, m_cbTypes);
  m_layout->addRow(m_buttons);
  
  // add the m_layout to the dialog
  this->setLayout(m_layout);
  
  // connect useful signals to slots
  connect(m_buttons, SIGNAL(accepted()), this, SLOT(btOkClicked()));
  connect(m_buttons, SIGNAL(rejected()), this, SLOT(btCancelClicked()));
  
  // the dialog is modal
  this->setModal(true);
  
  m_okClicked = false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TypeAndNameDialog::~TypeAndNameDialog()
{
  delete m_labConcreteType;
  delete m_labName;
  delete m_editName;
  delete m_cbTypes;
  delete m_layout;
  delete m_buttons;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


QString TypeAndNameDialog::show(const QStringList & types, QString & concreteType)
{
  // if the list is empty, there is no need to continue
  //  if(types.isEmpty())
  //   return QString();
  
  // clear the QComboBox and add the new m_items
  m_cbTypes->clear();
  m_cbTypes->addItems(types);
  
  // show the dialog (will not return while the dialog is visible)
  this->exec();
  
  // if the user did not clicked on "OK" or has not entered a name
  // then return an empty string
  if(!m_okClicked || m_editName->text().trimmed() == "")
    return QString();
  
  // set the selected concrete type and return the name
  concreteType = m_cbTypes->currentText();
  return m_editName->text();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TypeAndNameDialog::setName(const QString & newName)
{
  m_editName->setText(newName);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TypeAndNameDialog::btOkClicked()
{
  m_okClicked = true;
  this->setVisible(false); // this makes show() execution to continue
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TypeAndNameDialog::btCancelClicked()
{
  m_okClicked = false;
  this->setVisible(false); // this makes show() execution to continue
}
