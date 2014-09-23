#include <QtGui>

#include "ClientServer/client/CommitDetails.hh"

#include "ClientServer/client/CommitDetailsDialog.hh"

using namespace COOLFluiD::client;

CommitDetailsDialog::CommitDetailsDialog(QWidget * parent)
: QDialog(parent)
{
  //CommitDetails details;
  m_mainLayout = new QVBoxLayout(this);
  m_btOk = new QPushButton("Ok", this);
  m_buttonBox = new QDialogButtonBox(this);
  m_view = new QTableView(this);
  m_commitDetails = new CommitDetails();
  
  //m_view->setModel(&details);
  
  m_buttonBox->addButton(m_btOk, QDialogButtonBox::AcceptRole);
  
  m_view->resizeRowsToContents();
  //this->resize(this->width() * 2, this->height());
  
  m_mainLayout->addWidget(m_view, 0);
  m_mainLayout->addWidget(m_buttonBox, 1);
  
  this->setFixedSize(this->width(), this->height());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CommitDetailsDialog::~CommitDetailsDialog()
{
  
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetailsDialog::setCommitDetails(CommitDetails * details)
{
  m_commitDetails = details;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetailsDialog::show(CommitDetails & details)
{
  if (details.rowCount() > 0) 
  {
    m_view->setModel(&details);
    this->exec();  
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetailsDialog::show()
{
  this->show(*m_commitDetails);
}