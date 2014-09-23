#ifndef COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSvt_Ghost_hh
#define COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSvt_Ghost_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the NoSlipWall
   * Isothermal BC for coupled systems
   *
   * @author Thomas Wuilbaut
   *
   */

template <class MODEL>
class CoupledNoSlipWallIsothermalNSvt_Ghost : public NoSlipWallIsothermalNSvt<MODEL> {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledNoSlipWallIsothermalNSvt_Ghost(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CoupledNoSlipWallIsothermalNSvt_Ghost();

  /**
   * Configuration of the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /**
   * Create map between faceIdx and the data transfered
   */
  void setIndex();

private:

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet _sockets;

  //index of the data related to index of the face
  std::vector<CFint> _coupledDataID;

  // name of the datahandle containing the nodal values
  std::string _interfaceName;

  //back up value for the wall temperature imposed by user
  CFreal m_wallTempConst;

  bool _setIndex;

  //number of iterations during which to use the default values
  CFuint _defaultIterations;

}; // end of class CoupledNoSlipWallIsothermalNSvt_Ghost

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CoupledNoSlipWallIsothermalNSvt_Ghost.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSvt_Ghost_hh
