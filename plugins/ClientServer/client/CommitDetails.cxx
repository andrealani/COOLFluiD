#include <QtCore>

#include "ClientServer/client/CommitDetailsItem.hh"
#include "ClientServer/client/CommitDetails.hh"

using namespace COOLFluiD::client;

CommitDetails::CommitDetails(QObject * parent, const QString & nodePath)
: QAbstractItemModel(parent)
{
  m_nodePath = nodePath;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QVariant CommitDetails::data(const QModelIndex &index, int role) const
{
  QVariant returnValue; 
  
  
  if(role == Qt::DisplayRole)
  {
    int rowNumber = index.row();
    int colNumber = index.column();
    
    if(rowNumber >= 0 && rowNumber < m_items.size())
    {
      CommitDetailsItem * item = m_items.at(rowNumber);
      
      switch (colNumber) 
      {
        case 0:
          returnValue = item->getOptionName();
          break;
          
        case 1:
          returnValue = item->isNewOption() ? "Added" : "Modified";
          break;
          
        case 2:
        {
          QString oldValue = item->getOldValue();
          returnValue = oldValue.isEmpty() ? "--" : QString("\"%1\"").arg(oldValue);
          break;     
        }
          
        case 3:
        {
          QString value = item->getCurrentValue();
          returnValue = value.isEmpty() ? "--" : QString("\"%1\"").arg(value);
          break;     
        }
          
        default:
          break;
      }
      
    }   
  }
  
  return returnValue;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QVariant CommitDetails::headerData(int section, Qt::Orientation orientation,
                                   int role) const
{
  QVariant returnValue;
  
  if (role == Qt::DisplayRole) 
  {
    if(orientation == Qt::Horizontal)
    {
      switch (section) 
      {
        case 0:
          returnValue = "Name";
          break;
          
        case 1:
          returnValue = "Status";
          break;
          
        case 2:
          returnValue = "Old Value";
          break;     
          
        case 3:
          returnValue = "New value";
          break;
      }
    }
    else 
    {
      returnValue = QString("Option #%1").arg(section+1);
    }
  }
  return returnValue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex CommitDetails::index(int row, int column,
                                 const QModelIndex & parent) const
{
  CommitDetailsItem * item;
  QModelIndex index;
  
  if(!this->hasIndex(row, column, parent))
    return QModelIndex();
  
  if(m_items.isEmpty())
    return QModelIndex();
  
  item = m_items.at(row);
  index = createIndex(row, column, item);
  
  return index;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QModelIndex CommitDetails::parent(const QModelIndex &index) const
{
  if(!index.isValid())
    return QModelIndex();
  
  CommitDetailsItem * item = static_cast<CommitDetailsItem *> (index.internalPointer());
  
  if (item == NULL)
    return QModelIndex();
  
  return QModelIndex();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int CommitDetails::rowCount(const QModelIndex &parent) const
{
  return m_items.size();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int CommitDetails::columnCount(const QModelIndex &parent) const
{
  return 4;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetails::setOption(const QString & optionName, 
                              const QString & oldValue, 
                              const QString & currentValue)
{
  m_items << new CommitDetailsItem(optionName, oldValue, currentValue);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetails::setNewOption(const QString & optionName, 
                                 const QString & value, TOptionTypes type)
{ 
  m_items << new CommitDetailsItem(optionName, value);
} 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::setOptionNewValue(const QString & optionName, 
                                      const QString & value)
{
  if(!m_options.contains(optionName))
    return false;
  
  m_optionsNewValues[optionName] = value;
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::setOptionOldValue(const QString & optionName, 
                                      const QString & value)
{
  if(!m_options.contains(optionName))
    return false;
  
  m_optionsOldValues[optionName] = value;
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::setNewOptionValue(const QString & optionName, 
                                      const QString & value)
{
  if(!m_newOptions.contains(optionName))
    return false;
  
  m_newOptionsValues[optionName] = value;
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::removeOption(const QString & optionName)
{
  if(!m_options.contains(optionName))
    return false;
  
  m_options.removeOne(optionName);
  m_optionsOldValues.remove(optionName);
  m_optionsNewValues.remove(optionName);
  
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::removeNewOption(const QString & optionName)
{
  if(!m_newOptions.contains(optionName))
    return false;
  
  m_newOptions.removeOne(optionName);
  m_newOptionsValues.remove(optionName);
  m_newOptionsTypes.remove(optionName);
  
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::getOption(const QString & optionName, QString & oldValue, 
                              QString & newValue) const
{
  if(!m_options.contains(optionName))
    return false;
  
  oldValue = m_optionsOldValues[optionName];
  newValue = m_optionsNewValues[optionName];
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::getNewOption(const QString & optionName, 
                                 QString & value) const
{
  if(!m_newOptions.contains(optionName))
    return false;
  
  value = m_newOptionsValues[optionName];
  return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::contains(const QString & optionName, 
                             bool * isNewOption) const
{
  bool contains;
  
  contains = m_options.contains(optionName);
  
  if(contains && isNewOption != NULL)
    *isNewOption = false;
  
  else
  {
    contains = m_newOptions.contains(optionName);
    
    if(contains && isNewOption != NULL)
      *isNewOption = true;
  }
  
  return contains;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::isEmpty() const
{
  return m_options.isEmpty() && m_newOptions.isEmpty();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::hasOptions() const
{
  return m_options.isEmpty();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetails::hasNewOptions() const
{
  return m_newOptions.isEmpty();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetails::clear()
{
  this->clearOptions();
  this->clearNewOptions();
  m_nodePath.clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetails::clearOptions()
{
  m_options.clear();
  m_optionsNewValues.clear();
  m_optionsOldValues.clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetails::clearNewOptions()
{
  m_newOptions.clear();
  m_newOptionsValues.clear();
  m_newOptionsTypes.clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int CommitDetails::getOptionCount() const
{
  return m_options.count();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int CommitDetails::getNewOptionCount() const
{
  return m_newOptions.count();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TOptionTypes CommitDetails::getNewOptionType(const QString & optionName) const
{
  if(m_newOptions.contains(optionName))
    return m_newOptionsTypes[optionName];
  
  return NO_TYPE;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString CommitDetails::getNodePath() const
{
  return m_nodePath;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CommitDetails::setNodePath(const QString & nodePath)
{
  m_nodePath = nodePath;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString CommitDetails::toString() const
{
  QString details;
  QString oldValue;
  QString newValue;
  QString separator = QString("\n").rightJustified(30, '+');
  QString spaces = "   ";
  
  if(!m_options.isEmpty())
  {
    QStringList::const_iterator it = m_options.begin();
    
    details.append("Modified option(s):\n------------------------\n");
    
    while(it != m_options.end())
    {
      QString name = *it;
      
      oldValue = QString("Old value: \"%1\"\n").arg(m_optionsOldValues[name]);
      newValue = QString("New value: \"%1\"\n").arg(m_optionsNewValues[name]);
      
      details.append(QString("Option name: %1\n").arg(name));
      details.append(spaces + oldValue);
      details.append(spaces + newValue);
      
      it++;
      
      if(it != m_options.end())
        details.append(separator);
    }
  }  
  
  if(!m_newOptions.isEmpty())
  {
    if(!details.isEmpty())
      details.append(separator + separator);
    
    details.append("New option(s):\n-------------------\n");
    
    details.append(m_newOptions.join("\n" + spaces).prepend(spaces));
  }
  
  return details;
}
