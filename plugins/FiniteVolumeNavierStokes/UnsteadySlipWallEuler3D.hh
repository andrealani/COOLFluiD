#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadySlipWallEuler3D_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadySlipWallEuler3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the UnsteadySlipWall bc
   *
   * @author Thomas Wuilbaut
   *
   */
class UnsteadySlipWallEuler3D : public FVMCC_BC {

public:

  /**
   * Constructor
   */
  UnsteadySlipWallEuler3D(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadySlipWallEuler3D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private://function

  /**
   * Compute the mesh speed at the ghost state
   */
  void computeGhostStateSpeed(Framework::GeometricEntity *const face);

 private://data

  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// handle to the nodes at their past position
  Framework::DataSocketSink< Framework::Node*> socket_pastNodes;

  /// handle to the nodes at their future position
  Framework::DataSocketSink< Framework::Node*> socket_futureNodes;

  /// Temporary Storage of coordinates
  RealVector _bCoord;

  /// Speed of the mesh at the face
  RealVector _speed;

}; // end of class UnsteadySlipWallEuler3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadySlipWallEuler3D_hh
