#ifndef COOLFluiD_Numerics_FiniteVolume_CoupledPorousEuler2D_hh
#define COOLFluiD_Numerics_FiniteVolume_CoupledPorousEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/VectorialFunction.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DCons;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to implement a
 * coupled boundary condition for a porous wall
 * The data transfered should be in Euler2DCons variables
 *
 * @author Thomas Wuilbaut
 *
 */

class CoupledPorousEuler2D : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledPorousEuler2D(const std::string& name);

  /**
   * Default destructor
   */
  ~CoupledPorousEuler2D();

  /**
   * Configuration of the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  // physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DCons> _varSet;

  // physical model data
  RealVector _dataInnerState;

  // physical model data
  RealVector _dataGhostState;

  // name of the datahandle containing the nodal values
  std::string _interfaceName;

  // a vector of string to hold the functions
  std::vector<std::string> _functions;

  // a vector of string to hold the functions
  std::vector<std::string> _vars;

  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

}; // end of class CoupledPorousEuler2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CoupledPorousEuler2D_hh
