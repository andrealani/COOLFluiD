#ifndef COOLFluiD_Framework_SocketBundleSetter_hh
#define COOLFluiD_Framework_SocketBundleSetter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/TopologicalRegionSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Framework {

//////////////////////////////////////////////////////////////////////////////

class SocketBundle {
public:
  
  SocketBundle();
  
  /// the socket to the radiative heat flux at the wall faces
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > states;
  
  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> gstates;
  
  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > nodes;
  
  /// storage of the nodal state's
  Framework::DataSocketSink < RealVector > nstates;
  
  /// handle to the face normals
  Framework::DataSocketSink< CFreal> normals;
  
  /// IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSink<CFint> isOutward;
  
  /// IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSink<CFreal> volumes;
  
  /// storage of face centroids
  Framework::DataSocketSink<CFreal> faceCenters;
  
  /// storage of face areas
  Framework::DataSocketSink<CFreal> faceAreas; 
  
  /// storage of the binned opacity
  Framework::DataHandle<CFreal> alpha_avbin;
  
  /// storage of the binned radiative source
  Framework::DataHandle<CFreal> B_bin;
  
};

//////////////////////////////////////////////////////////////////////////////

class SocketBundleSetter {
public:

  SocketBundleSetter();
  
  /// Default destructor
  virtual ~SocketBundleSetter();
  
  /// set up private data
  void setDataSockets(SocketBundle sockets);

  SocketBundle* getDataSocket(){return &m_sockets;}

  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder>*
       getCellTRSbuilder(){return &m_cellBuilder;}

  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder>*
       getFaceTrsBuilder(){return &m_faceBuilder;}

protected:
  SocketBundle m_sockets;

  /// cell builder
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_cellBuilder;

  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceBuilder;

}; // end of class

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SocketBundleSetter_hh
