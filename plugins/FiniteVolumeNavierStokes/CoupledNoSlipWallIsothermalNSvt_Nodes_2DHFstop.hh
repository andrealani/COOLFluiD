#ifndef COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF_hh
#define COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

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
class CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF : public NoSlipWallIsothermalNSvt<MODEL> {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF();

  /**
   * Configuration of the command
   */
  virtual void configure (const Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Deallocate private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void unsetup();

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

  /// index of the data related to index of the node
  std::vector<CFint> _coupledDataID;

  /// Map to get back the index of the node in the TRS list from its LocalID
  Common::CFMap<CFuint, CFuint> _trsNodeIDMap;

  /// name of the datahandle containing the nodal values
  std::string _interfaceName;

  /// back up value for the wall temperature imposed by user
  CFreal m_wallTempConst;

  /// ???
  bool _setIndex;

  /// number of iterations during which to use the default values
  CFuint _defaultIterations;

  /// reference temp for adimensionalization
  CFreal m_refTemp;

protected:

  /// acquaintance of the concrete variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokes2DVarSet> _diffusiveVarSet;

  /// dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;

  /// storage of norm of the delta flux across interface
  CFreal m_sqr_delta_flux;


}; // end of class CoupledNoSlipWallIsothermalNSvt_Nodes

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CoupledNoSlipWallIsothermalNSvt_Nodes_2DHF_hh
