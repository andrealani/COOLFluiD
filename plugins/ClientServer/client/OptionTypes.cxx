#include <QtCore>

#include "ClientServer/client/OptionTypes.hh"

using namespace COOLFluiD::client;

QHash<TOptionTypes, QString> OptionTypes::types;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool OptionTypes::isValid(TOptionTypes id)
{
  OptionTypes::buildTypes();
  
  return OptionTypes::types.contains(id);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TOptionTypes OptionTypes::getTypeId(const QString & type)
{
  QHash<TOptionTypes, QString>::iterator it;
  
  OptionTypes::buildTypes();
  
  it = OptionTypes::types.begin();
  
  while(it != OptionTypes::types.end())
  {
    if(it.value() == type)
      return it.key();
    it++;
  }
  
  return NO_TYPE;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString OptionTypes::getTypeString(TOptionTypes type)
{
  OptionTypes::buildTypes();
  
  if(OptionTypes::isValid(type))
    return OptionTypes::types[type];
  
  return QString();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void OptionTypes::buildTypes()
{
  static bool mapBuilt = false;
  
  if(mapBuilt) // if the map has already been built...
    return;     // the function returns (there no need to build it again)
  
  OptionTypes::types[ TYPE_BOOL ] = "bool";
  OptionTypes::types[ TYPE_INT ] = "int";
  OptionTypes::types[ TYPE_UNSIGNED_INT ] = "unsigned int";
  OptionTypes::types[ TYPE_DOUBLE ] = "double";
  OptionTypes::types[ TYPE_STRING ] = "std::string";
  OptionTypes::types[ TYPE_FILES ] = "files";
  OptionTypes::types[ TYPE_LIBRARIES ] = "libraries";
  OptionTypes::types[ TYPE_HOST_LIST ] = "hostList";
  
  mapBuilt = true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QStringList OptionTypes::getTypesList()
{
  static QStringList list;
  static bool listBuilt = false;
  QHash<TOptionTypes, QString>::iterator it;
  
  if(listBuilt)   // if the list has already been built...
    return list;   // the function returns (there no need to build it again)
  
  OptionTypes::buildTypes();
  
  it = OptionTypes::types.begin();
  
  while(it != OptionTypes::types.end())
  {
    list << it.value();
    it++;
  }
  
  listBuilt = true;
  
  return list;
}
