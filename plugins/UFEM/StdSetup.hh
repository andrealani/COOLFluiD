#ifndef COOLFluiD_UFEM_StdSetup_hh
#define COOLFluiD_UFEM_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMSolverData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to setup the MeshData.
/// @author Tiago Quintino
class UFEM_API StdSetup : public UFEMSolverCom {
public:

  /// Constructor.
  explicit StdSetup(const std::string& name);

  /// Destructor.
  virtual ~StdSetup();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Execute Processing actions
  void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// Set the mapping between each boundary Face and its neighbor
  /// cell, updated version
  void setBndFace2InnerCellAndBndFace2BndStateConnectivity();

  /// Build info about where the closest wall is located.
  /// This is just a small thing that reads in the data
  /// from a file. It is intended to use for debug purposes only.
  void setWallNearestSegment();
  
  /// Read Wall distance from file
  void readWallDistanceFromFile();

protected:

  /// mapping nodeID to stateID
  std::vector<CFuint> _nodeIdToStateId;

  /// socket with proxy to be able to use the data handle of nodal states uniformly
  /// independently from the actual storage type being RealVector or State*
  Framework::DataSocketSource<Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;

  /// socket with flags to check if a state's equation has been updated or not
  Framework::DataSocketSource<bool> socket_isUpdated;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// socket for GeometricEntity mapper to corresponding TopologicalRegionSet 
  /// (boundary trs geoentities and wall adjacent inner cell geoentities connectivity)
  Framework::DataSocketSource< std::vector<CFuint> > socket_connBndFace2InnerCell;

  /// socket for boundary trs-local connectivity table between local and trs-local numbering of states
  /// valid only for non-innercell trs
  Framework::DataSocketSource< Common::ConnectivityTable<CFuint> > socket_connBndFace2BndState;

  /// which trs is the wall
  CFuint wallNearestTrs;

  /// index of the closest wall boundary segment
  Framework::DataSocketSource< CFuint > socket_wallNearestSegment;

  /// distance to the closest wall boundary segment
  Framework::DataSocketSource< CFreal > socket_wallNearestDistance;

  /// distance to the closest wall boundary segment, nodal wise
  Framework::DataSocketSource< CFreal > socket_wallNearestDistanceState;

  /// wall-velocity gradient
  Framework::DataSocketSource< CFreal > socket_wallNearestVelocityGradient;

  /// the socket to the data handle of arrays of flags specifying if a strong BC has been applied
  /// for the given variables in boundary states
  Framework::DataSocketSource<std::vector<bool> > socket_appliedStrongBC;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_StdSetup_hh

