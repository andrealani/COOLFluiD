#ifndef COOLFluiD_Numerics_MeshLaplacianSmoothing_StdPrepare_hh
#define COOLFluiD_Numerics_MeshLaplacianSmoothing_StdPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "LaplacianSmoothingData.hh"
#include "MeshTools/QualityCalculator.hh"
#include "Common/SafePtr.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */

class StdPrepare : public LaplacianSmoothingCom {
public:

  /**
   * Constructor.
   */
  StdPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~StdPrepare()
  {
  }

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

  /**
   * For each node, compute the quality of the worse cell connected to it
   */
  void computeWorstQualityCells();

private: // data

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;

  /// storage of displacements vector
  Framework::DataSocketSink<RealVector> socket_nodalDisplacements;

  /// the socket to the data handle of the node connectivity
  Framework::DataSocketSink<std::vector<CFuint> > socket_nodeToCellConnectivity;

  /// storage of "quality at a node"
  Framework::DataSocketSink<CFreal> socket_qualityNode;

  /// handle for the cells storage
  Common::SafePtr<std::vector<Framework::GeometricEntity*> > _cells;

  /// Functor that computes the requested quality
  Common::SelfRegistPtr<MeshTools::QualityCalculator> _computeQuality;

}; // class StdPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshLaplacianSmoothing_StdPrepare_hh

