#ifndef COOLFluiD_Numerics_MeshFEMMove_OscillatingAirfoilPrepare_hh
#define COOLFluiD_Numerics_MeshFEMMove_OscillatingAirfoilPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "BasePrepare.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to define
   * the movement of the boundaries
   *
   * @author Thomas Wuilbaut
   *
   */

class OscillatingAirfoilPrepare : public BasePrepare {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  OscillatingAirfoilPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~OscillatingAirfoilPrepare()
  {
  }

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up the data
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private: // functions

  /**
   * Move the boundary nodes
   */
  virtual void moveBoundaries();

private: // data

  // the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  // Angle of the OscillatingAirfoil
  CFreal _currentAlpha;

  // rotation center for config
  std::vector<CFreal> _rotationCenterConf;

  // rotation center
  RealVector _rotationCenter;

}; // class OscillatingAirfoilPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_OscillatingAirfoilPrepare_hh

