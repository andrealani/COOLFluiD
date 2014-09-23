#include <QtCore>

#include "ClientServer/client/CommitDetailsItem.hh"

using namespace COOLFluiD::client;

CommitDetailsItem::CommitDetailsItem(const QString & optionName, 
                                     const QString & oldValue, 
                                     const QString & currentValue)
{
  m_optionName = optionName;
  m_oldValue = oldValue;
  m_currentValue = currentValue;
  m_newOption = false;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CommitDetailsItem::CommitDetailsItem(const QString & oldValue, 
                                     const QString & currentValue)
{
  m_optionName = m_optionName;
  m_oldValue = oldValue;
  m_newOption = true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CommitDetailsItem::isNewOption() const
{
  return m_newOption;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString CommitDetailsItem::getCurrentValue() const
{
  return m_currentValue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString CommitDetailsItem::getOldValue() const
{
  return m_oldValue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QString CommitDetailsItem::getOptionName() const
{
  return m_optionName;
}
