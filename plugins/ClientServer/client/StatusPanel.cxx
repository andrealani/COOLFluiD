#include <QtGui>

#include "ClientServer/client/StatusModel.hh"

#include "ClientServer/client/StatusPanel.hh"

using namespace COOLFluiD::client;

StatusPanel::StatusPanel(StatusModel * model, QWidget * parent)
: QTreeView(parent)
{
  m_model = model;
  this->setModel(m_model);
  this->header()->setResizeMode(QHeaderView::ResizeToContents);
  
  connect(m_model, SIGNAL(subSystemAdded(const QModelIndex &)),
          this, SLOT(subSystemAdded(const QModelIndex &)));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

StatusPanel::~StatusPanel()
{
  delete m_model;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void StatusPanel::subSystemAdded(const QModelIndex & index)
{
  this->expand(index.parent());
  this->expand(index);
}
