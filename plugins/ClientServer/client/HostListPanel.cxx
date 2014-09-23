#include <QtCore>
#include <QtGui>
#include <QtXml>

#include "ClientServer/client/HostListPanel.hh"

using namespace COOLFluiD::client;

HostListPanel::HostListPanel()
: QListView()
{
  m_viewModel = new QStandardItemModel();
  
  this->setModel(m_viewModel);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

HostListPanel::~HostListPanel()
{
  this->clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void HostListPanel::setItems(const QDomDocument & doc)
{
  QDomNodeList childNodes = doc.firstChild().childNodes();
  
  this->clear();
  
  for(int i = 0 ; i < childNodes.count() ; i++)
  {
    QDomNode node = childNodes.at(i);
    QDomNamedNodeMap attrs = node.attributes();
    QDomText text = node.firstChild().toText();
    QStandardItem * item = new QStandardItem(text.nodeValue());
    bool selected = node.attributes().namedItem("selected").nodeValue() == "true";
    
    item->setCheckable(true);
    item->setEditable(false);
    
    item->setCheckState(selected ? Qt::Checked : Qt::Unchecked);
    
    m_viewItems.append(item);
    m_viewModel->appendRow(item);
  }
  
  m_hosts = doc;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QDomDocument HostListPanel::getDocument() const
{
  QDomDocument doc = m_hosts.cloneNode(true).toDocument();
  QDomNodeList childNodes = doc.firstChild().childNodes();
  
  for(int i = 0 ; i < childNodes.size() ; i++)
  {
    QDomElement node = childNodes.at(i).toElement();
    bool selected = m_viewItems.at(i)->checkState() == Qt::Checked;
    
    node.setAttribute("selected", selected ? "true" : "false");
  }
  
  return doc;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool HostListPanel::isModified() const
{
  QDomNodeList childNodes = m_hosts.firstChild().childNodes();
  bool modified = false;
  
  for(int i = 0 ; i < childNodes.size() && !modified; i++)
  {
    QDomElement node = childNodes.at(i).toElement();
    QString selected;
    
    if(m_viewItems.at(i)->checkState() == Qt::Checked)
      selected = "true";
    else
      selected = "false";
    
    modified = node.attribute("selected") != selected;
    
  }
  
  return modified;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void HostListPanel::clear()
{
  m_hosts.clear();
  m_viewModel->clear();
  
  // list m_items were deleted (memory freed) when clearing the model, 
  // so we just have to clear the list
  m_viewItems.clear();
}