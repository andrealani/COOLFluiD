#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_UpdateMesh_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_UpdateMesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpringAnalogyData.hh"
#include "Framework/Storage.hh"
#include "Framework/Node.hh"
#include "BallVertexCalculator.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order adapt the Mesh.
   *
   * @author Thomas Wuilbaut
   *
   */

class UpdateMesh : public SpringAnalogyCom {
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

  /**
   * Move the node
   */
  void computeNewCoordinates(const CFuint iNode);

  /**
   * Compute the vector for displacement of the node
   */
  void computeDisplacementVectors();

  /**
   * Check the Volume of the newly created cells
   */
  void checkNewCells();

  // Compute the angle at the node
  //CFreal computeNodalAngle(std::vector<Framework::Node*>* cellNodes, CFuint iNode);

private: //data

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  ///storage for the average vector
  Framework::DataSocketSink<RealVector> socket_averageVector;

  /// storage of the sum of the weighting factors
  Framework::DataSocketSink<CFreal> socket_sumWeight;

  /// storage that indicates if nodes are movable
  Framework::DataSocketSink<bool> socket_isMovable;

  /// storage of displacements vector
  Framework::DataSocketSink<RealVector> socket_nodalDisplacements;

  /// storage of value of the distance from the node to the closest node of the boundary
  Framework::DataSocketSink<CFreal> socket_wallDistance;

  ///temporary RealVector
  RealVector _coordI;

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

  CFuint _iIter;

  //the computer of the ball vertex spring
  BallVertexCalculator _ballVertexComputer;

}; // class UpdateMesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_UpdateMesh_hh

