
#ifndef COOLFluiD_Numerics_FluctSplit_VorticitySetup_hh
#define COOLFluiD_Numerics_FluctSplit_VorticitySetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// to be executed in order to create data for the FluctSplit Method
/// @author Tiago Quintino
/// @author Andrea Lani
class FluctSplit_API VorticitySetup : public FluctuationSplitCom {
public: // member functions

  /// Constructor.
  explicit VorticitySetup(const std::string& name);

  /// Destructor.
  virtual ~VorticitySetup();

  /// Configure the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

//   /// Execute Processing actions
  virtual void execute();

/// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
   static void defineConfigOptions(Config::OptionList& options);
private: // helper member functions

  /// Compute the boundary data
  void computeNormalsData();

  /// Flag the boundary states
  void flagBoundaryStates();

  /// Set the mapping between each boundary Face and its neighbor
  /// cell
  void setFaceNeighCell();

  /// Store the (outward) normals at the boundaries
  void storeBoundaryNormals();

  /// Compute the cell volumes
  void computeCellVolumes();

protected: // data

  /// mapping from nodeID to stateID
  std::vector<CFuint> m_nodeIdToStateId;

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_dynamicSockets;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for inward normals
  Framework::DataSocketSource< InwardNormalsData*> socket_normals;

  /// socket with proxy to be able to use the data handle of nodal states uniformly
  /// independently from the actual storage type being RealVector or State*
  Framework::DataSocketSource
    <Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;

  /// socket with flags to check if a state has been updated or not
  Framework::DataSocketSource<bool> socket_isUpdated;

  /// socket with flags to check if a state is on the boundary
  Framework::DataSocketSource<bool> socket_isBState;

  // the socket to the data handle of arrays of flags specifying if the
  // time jacobian contribution of certain variables in boundary states
  // have to be discarded
  Framework::DataSocketSource<std::vector<bool> > socket_discardTimeJacob;

  /// socket for inward normals data
  Framework::DataSocketSource<CFreal> socket_normalsData;

  /// handle to the neighbor cell
  Framework::DataSocketSource<
  Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// socket with tempSize
  Framework::DataSocketSource<CFuint> socket_tempSize;

  /// socket with volumes
  Framework::DataSocketSource<CFreal> socket_volumes;

//   /// socket with theta clending coefficient
       Framework::DataSocketSource<CFreal> socket_vorticity;
  // /// initial limiter socket name
   std::string _vorticitySocketName;
   


  /// storage for flag
  Framework::DataSocketSource<CFuint> socket_flagStates;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdSetup_hh

