#include <QString>

#include "ClientServer/client/GlobalLog.hh"

using namespace COOLFluiD::client;

GlobalLog * GlobalLog::m_instance = NULL;

GlobalLog::GlobalLog()
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GlobalLog * GlobalLog::getInstance()
{
  if(m_instance == NULL)
    m_instance = new GlobalLog();
  
  return m_instance;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GlobalLog::error(const QString & msg, bool fromServer)
{
  getInstance()->writeError(msg, fromServer);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GlobalLog::message(const QString & msg, bool fromServer)
{
  getInstance()->writeMessage(msg, fromServer);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GlobalLog::writeError(const QString & msg, bool fromServer) const
{
  if(!msg.isEmpty())
    emit sigError(msg, fromServer);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GlobalLog::writeMessage(const QString & msg, bool fromServer) const
{
  if(!msg.isEmpty())
    emit sigMessage(msg, fromServer);
}
