#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSPvtSimple_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSPvtSimple_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
    
    namespace NavierStokes {
      class EulerVarSet;
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
class NoSlipWallIsothermalNSPvtSimple : public FVMCC_BC { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  NoSlipWallIsothermalNSPvtSimple(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNSPvtSimple();
  
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
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _updateVarSet;
  
  /// temperature at the wall
  CFreal _wallTemp;
    
}; // end of class NoSlipWallIsothermalNSPvtSimple

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSPvtSimple_hh
