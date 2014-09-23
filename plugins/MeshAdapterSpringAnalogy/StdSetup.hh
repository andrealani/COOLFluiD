#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_StdSetup_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpringAnalogyData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to setup the MeshData.
 *
 * @author Thomas Wuilbaut
 *
 */
class StdSetup : public SpringAnalogyCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(std::string name);

  /**
   * Destructor.
   */
  ~StdSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private: // functions

  /// Create the list of Nodes of the TRSs
  void createNodesListInTRS();

  /// Flag the boundary nodes
  void flagBoundaryNodes();

  /// Create the Node to Cell Connectivity
  void createNodeToCellConnectivity();

protected: // data

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;


  ///storage for the average vector
  Framework::DataSocketSource<RealVector> socket_averageVector;

  /// storage of the sum of the weighting factors
  Framework::DataSocketSource<CFreal> socket_sumWeight;

  /// storage that indicates if nodes are movable
  Framework::DataSocketSource<bool> socket_isMovable;

  /// storage of displacements vector
  Framework::DataSocketSource<RealVector> socket_nodalDisplacements;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_StdSetup_hh

