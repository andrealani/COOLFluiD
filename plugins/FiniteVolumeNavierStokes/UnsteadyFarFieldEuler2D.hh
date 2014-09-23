#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadyFarFieldEuler2D_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadyFarFieldEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to implement UnsteadyFarField boundary condition
 *
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 *
 */
class UnsteadyFarFieldEuler2D : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  UnsteadyFarFieldEuler2D(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadyFarFieldEuler2D();

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

private: // function


  /**
   * Compute the mesh speed at the ghost state
   */
  void computeBoundarySpeed(Framework::GeometricEntity *const face);

 private: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// handle to the nodes at their past position
  Framework::DataSocketSink< Framework::Node*> socket_pastNodes;

  /// handle to the nodes at their future position
  Framework::DataSocketSink< Framework::Node*> socket_futureNodes;

  /// Temporary Storage of coordinates
  RealVector _coord;

  /// Speed of the mesh at the face
  RealVector _speed;

  /// static temperature
  CFreal _temperature;

  /// static pressure
  CFreal _pressure;

  /// x velocity
  CFreal _uInf;

  /// y velocity
  CFreal _vInf;

}; // end of class UnsteadyFarFieldEuler2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadyFarFieldEuler2D_hh
