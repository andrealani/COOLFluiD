#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNS2D_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNS2D_hh

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
   * This class represents a command that applies the mirror bc 
   * 
   * @author Andrea Lani 
   *
   */
class NoSlipWallAdiabaticNS2D : public FVMCC_BC { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  NoSlipWallAdiabaticNS2D(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NoSlipWallAdiabaticNS2D();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private:
  
  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;
  
  /// physical model data
  RealVector _dataInnerState;
  
  /// physical model data
  RealVector _dataGhostState;
  
  /// X-component of a velocity vector of the wall
  CFreal _xWallVelocity;
    
  /// Y-component of a velocity vector of the wall
  CFreal _yWallVelocity;
    
}; // end of class NoSlipWallAdiabaticNS2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNS2D_hh
