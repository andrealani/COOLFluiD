#ifndef COOLFluiD_Numerics_MeshLaplacianSmoothing_StdSetup_hh
#define COOLFluiD_Numerics_MeshLaplacianSmoothing_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "LaplacianSmoothingData.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand
 * action to be executed in order to setup the LaplacianSmoothing Method
 *
 * @author Thomas Wuilbaut
 */
class StdSetup : public LaplacianSmoothingCom {
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

private: // data

  /// Socket for the average vector
  Framework::DataSocketSource<RealVector> socket_averageVector;

  /// Socket for the sum of the weighting factors
  Framework::DataSocketSource<CFreal> socket_sumWeight;

  /// Socket with values of the quality associated to a node.
  /// Quality of the worse cell connected to the node.
  Framework::DataSocketSource<CFreal> socket_qualityNode;

  /// storage of displacements vector
  Framework::DataSocketSource<RealVector> socket_nodalDisplacements;

  /// the socket to the data handle of the node connectivity
  Framework::DataSocketSource<std::vector<CFuint> > socket_nodeToCellConnectivity;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshLaplacianSmoothing_StdSetup_hh

