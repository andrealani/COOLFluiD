#include <QtCore>
#include <stdexcept>

#include "ClientServer/client/UnknownTypeException.hh"

#include "ClientServer/client/FilesListItem.hh"

using namespace COOLFluiD::client;

FilesListItem::FilesListItem(const QIcon & icon, const QString & text,
                             FilesListItemType type)
: QStandardItem(icon, text)
{
  if(type != DIRECTORY && type != FILE)
    throw UnknownTypeException(FromHere(), "Unknown item type");
  
  m_type = type;
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FilesListItemType FilesListItem::getType() const
{
  return m_type;
}
