#ifndef COOLFluiD_Numerics_MeshLaplacianSmoothing_UpdateMesh_hh
#define COOLFluiD_Numerics_MeshLaplacianSmoothing_UpdateMesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "LaplacianSmoothingData.hh"
#include "Framework/Storage.hh"
#include "Framework/Node.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order adapt the Mesh.
   *
   * @author Thomas Wuilbaut
   *
   */

class UpdateMesh : public LaplacianSmoothingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  UpdateMesh(const std::string& name);

  /**
   * Destructor.
   */
  ~UpdateMesh()
  {
  }

  /**
   * set Up member data
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private: //functions

  void computeNewCoordinates(const CFuint iNode);

  /**
   * Check the Volume of the newly created cells
   */
  void checkNewCells();

  /**
   * Compute the average displacement vector
   */
  void computeDisplacementVectors();

private: //data

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  ///storage for the average vector
  Framework::DataSocketSink<RealVector> socket_averageVector;

  /// storage of the sum of the weighting factors
  Framework::DataSocketSink<CFreal> socket_sumWeight;

  /// storage of value of the distance from the node to the closest node of the boundary
  Framework::DataSocketSink<CFreal> socket_wallDistance;

  /// storage of the nodal quality
  Framework::DataSocketSink<CFreal> socket_qualityNode;

  /// storage of the node to cell connectivity
  Framework::DataSocketSink<std::vector<CFuint> > socket_nodeToCellConnectivity;

  /// Relaxation factor
  CFreal _relaxation;

  /// Number of Smoothing Iterations
  CFuint _nbSmoothingIter;

  /// Type of weighting function to use
  std::string _weightType;

  /// threshold for selection of nodes/cells to modify
  CFreal _qualityThreshold;

  // Flag for at least one cell with negative Volume
  bool _hasNegativeVolumeCells;

  // Number of Cells with a negative volume
  CFuint _nbNegativeVolumeCells;

  //temporary RealVector
  RealVector _coordI;

  CFuint _iIter;
}; // class UpdateMesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshLaplacianSmoothing_UpdateMesh_hh

