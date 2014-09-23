#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadySubInletEuler2DUVT_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadySubInletEuler2DUVT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/Node.hh"
#include "Framework/DataSocketSink.hh"

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
   * This class represents a subsonic inlet command with the
   * initial conditions given for tTotal, pTotal and alpha
   *
   * @author Mehmet Sarp Yalim
   * @author Andrea Lani
   * @author Thomas Wuilbaut
   */
class UnsteadySubInletEuler2DUVT : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  UnsteadySubInletEuler2DUVT(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadySubInletEuler2DUVT();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private://function

  /**
   * Compute the mesh speed at the ghost state
   */
  void computeGhostStateSpeed(Framework::GeometricEntity *const face);

 private: // data

  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// x velocity
  CFreal                                   _uinf;

  /// y velocity
  CFreal                                   _vinf;

  /// static temperature
  CFreal                                   _temperature;

  /// handle to the nodes at their past position
  Framework::DataSocketSink< Framework::Node*> socket_pastNodes;

  /// handle to the nodes at their future position
  Framework::DataSocketSink< Framework::Node*> socket_futureNodes;

  /// Storage of coordinates
  RealVector _bCoord;

  /// Speed of the mesh at the face
  RealVector _speed;

  /// Vector for U,V and T computed from the function
  RealVector _dimState;

  /// Vector for coordinates + time
  RealVector _variables;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

}; // end of class UnsteadySubInletEuler2DUVT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadySubInletEuler2DUVT_hh
