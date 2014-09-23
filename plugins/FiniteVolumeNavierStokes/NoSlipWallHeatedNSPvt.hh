#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallHeatedNSPvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallHeatedNSPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies a heated wall bc
   *
   * @author Thomas Wuilbaut
   *
   */
class NoSlipWallHeatedNSPvt : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallHeatedNSPvt(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallHeatedNSPvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffusiveVarSet;
  
  /// X-component of a velocity vector of the wall
  CFreal _xWallVelocity;

  /// Y-component of a velocity vector of the wall
  CFreal _yWallVelocity;

  /// Z-component of a velocity vector of the wall
  CFreal _zWallVelocity;
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;
  
  ///Normal Heat Flux to impose
  CFreal _heatFlux;
  
}; // end of class NoSlipWallHeatedNSPvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallHeatedNSPvt_hh
