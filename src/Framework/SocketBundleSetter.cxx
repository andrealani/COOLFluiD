#include "Framework/SocketBundleSetter.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Framework {

//////////////////////////////////////////////////////////////////////////////

SocketBundle::SocketBundle():
  states("Null"),
  gstates("Null"),
  nodes("Null"),
  nstates("Null"),
  normals("Null"),
  isOutward("Null"),
  volumes("Null"),
  faceCenters("Null"),
  faceAreas("Null"),
  alpha_avbin(CFNULL),
  B_bin(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

SocketBundleSetter::SocketBundleSetter()
{
}

//////////////////////////////////////////////////////////////////////////////

SocketBundleSetter::~SocketBundleSetter()
{
}

//////////////////////////////////////////////////////////////////////////////

void SocketBundleSetter::setDataSockets(SocketBundle sockets)
{
  m_sockets = sockets;
  // cell builder initialization
  m_cellBuilder.setup();
  m_cellBuilder.getGeoBuilder()->setDataSockets(sockets.states, sockets.gstates, sockets.nodes);
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  cellData.trs = MeshDataStack::getActive()->getTrs("InnerCells");

  // face builder initializaxtion
  m_faceBuilder.setup();
  m_faceBuilder.getGeoBuilder()->setDataSockets(sockets.states, sockets.gstates, sockets.nodes);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
