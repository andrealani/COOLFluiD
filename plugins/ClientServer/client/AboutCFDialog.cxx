#include <QtGui>

#include "Common/COOLFluiD.hh"

#include "ClientServer/client/AboutCFDialog.hh"

using namespace COOLFluiD::client;

AboutCFDialog::AboutCFDialog(QWidget * parent)
{
  this->setWindowTitle("About COOLFluiD");
  
  m_mainLayout = new QVBoxLayout(this);
  m_infoLayout = new QFormLayout();
  
  m_btOK = new QPushButton("OK");
  
  m_infoList << new CFInfo("COOLFluiD version:", COOLFLUID_VERSION_STR, m_infoLayout);
  m_infoList << new CFInfo("Kernel version:", CF_KERNEL_VERSION_STR, m_infoLayout);
  m_infoList << new CFInfo("Build operating system:", QString("%1 [%2bits]").arg(CF_OS_LONGNAME).arg(CF_OS_BITS), m_infoLayout);
  m_infoList << new CFInfo("Build processor:", CF_BUILD_PROCESSOR, m_infoLayout);
  
  m_mainLayout->addLayout(m_infoLayout);
  m_mainLayout->addWidget(m_btOK);
  
  connect(m_btOK, SIGNAL(clicked()), this, SLOT(accept()));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AboutCFDialog::~AboutCFDialog()
{
  //  delete m_btOK;
  //  delete m_mainLayout;
  //  delete m_infoLayout;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AboutCFDialog::CFInfo::CFInfo(const QString & name, const QString & value, QFormLayout * parent)
{
  this->labName = new QLabel(name);
  this->labValue = new QLabel(value);
  
  parent->addRow(this->labName, this->labValue);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AboutCFDialog::CFInfo::~CFInfo()
{
  delete this->labName;
  delete this->labValue;
}
