#ifndef COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSrvt_Nodes_hh
#define COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSrvt_Nodes_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSrvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the Isothermal wall
   * BC in case of update variables including velocity components and temperature
   * [p v T y_i] or [rho_i v T] etc.
   *
   * @author Andrea Lani
   *
   */
template <class MODEL>
class CoupledNoSlipWallIsothermalNSrvt_Nodes : public NoSlipWallIsothermalNSrvt<MODEL> {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledNoSlipWallIsothermalNSrvt_Nodes(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CoupledNoSlipWallIsothermalNSrvt_Nodes();

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
  Framework::DynamicDataSocketSet<> _sockets;

  //index of the data related to index of the node
  std::vector<CFint> _coupledDataID;

  // Map to get back the index of the node in the TRS list from its LocalID
  Common::CFMap<CFuint, CFuint> _trsNodeIDMap;

  // name of the datahandle containing the nodal values
  std::string _interfaceName;

  //back up value for the wall temperature imposed by user
  CFreal m_wallTempConst;

  bool _setIndex;

  //number of iterations during which to use the default values
  CFuint _defaultIterations;

  //reference temp for adimensionalization
  CFreal m_refTemp;

}; // end of class NoSlipWallIsothermalNSPvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CoupledNoSlipWallIsothermalNSrvt_Nodes.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSrvt_Nodes_hh
