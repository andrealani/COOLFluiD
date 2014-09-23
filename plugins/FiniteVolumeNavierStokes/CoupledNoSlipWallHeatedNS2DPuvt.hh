#ifndef COOLFluiD_Numerics_FiniteVolumeNavierStokes_CoupledNoSlipWallHeatedNS2DPuvt_hh
#define COOLFluiD_Numerics_FiniteVolumeNavierStokes_CoupledNoSlipWallHeatedNS2DPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallHeatedNSPvt.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies a heated wall bc
   * in a coupled simulation
   *
   * @author Thomas Wuilbaut
   *
   */
class CoupledNoSlipWallHeatedNS2DPuvt : public NoSlipWallHeatedNSPvt {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledNoSlipWallHeatedNS2DPuvt(const std::string& name);

  /**
   * Default destructor
   */
  ~CoupledNoSlipWallHeatedNS2DPuvt();

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
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

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

  //check if function set Index has been done
  bool _setIndex;

  //number of iterations during which to use the default values
  CFuint _defaultIterations;

}; // end of class CoupledNoSlipWallHeatedNS2DPuvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeNavierStokes_CoupledNoSlipWallHeatedNS2DPuvt_hh
