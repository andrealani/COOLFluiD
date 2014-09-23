#ifndef COOLFluiD_Numerics_FiniteElement_StdSetup_hh
#define COOLFluiD_Numerics_FiniteElement_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/TopologicalRegionSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to setup the MeshData.
 *
 * @author Tiago Quintino
 */
class StdSetup : public FiniteElementMethodCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /**
   * Set the mapping between each boundary Face and its neighbor
   * cell
   */
  void setFaceNeighCell();

protected:

  /// mapping nodeID to stateID
  std::vector<CFuint> _nodeIdToStateId;

  /// socket with proxy to be able to use the data handle of nodal states uniformly
  /// independently from the actual storage type being RealVector or State*
  Framework::DataSocketSource
    <Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;

  /// socket with flags to check if a state has been updated or not
  Framework::DataSocketSource<bool> socket_isUpdated;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  // the socket to the data handle of the node's
  Framework::DataSocketSink <Framework::Node* , Framework::GLOBAL> socket_nodes;
  
  /// socket for GeometricEntity mapper to corresponding TopologicalRegionSet
  Framework::DataSocketSource<
			      Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
                              socket_faceNeighCell;

  // the socket to the data handle of arrays of flags specifying if a strong BC has been applied
  // for the given variables in boundary states
  Framework::DataSocketSource<std::vector<bool> > socket_appliedStrongBC;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_StdSetup_hh

