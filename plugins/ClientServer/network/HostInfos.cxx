#include <QString>

#include "ClientServer/network/HostInfos.hh"

using namespace COOLFluiD::network;

HostInfos::HostInfos(const QString & hostname, int nbSlots, int maxSlots)
{
  m_hostname = hostname;
  m_nbSlots = nbSlots;
  m_maxSlots = maxSlots;
}
