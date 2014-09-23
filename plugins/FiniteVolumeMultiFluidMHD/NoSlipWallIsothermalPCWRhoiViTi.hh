#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalPCWRhoiViTi_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalPCWRhoiViTi_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MultiFluidMHD {      
      class DiffMFMHD2DVarSet;
    }
    
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a No Slip condition imposing the temperatures at the wall. 
   * Perfectly conducting condition for Maxwell equations
   * 
   * @author Alejandro Alvarez
   *
   */
class NoSlipWallIsothermalPCWRhoiViTi : public FVMCC_BC { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  NoSlipWallIsothermalPCWRhoiViTi(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalPCWRhoiViTi();
  
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
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > _updateVarSet;
  
  /// temperatures at the wall
  std::vector<CFreal> _wallTemp;
  
    
}; // end of class NoSlipWallIsothermalPCWRhoiViTi

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalPCWRhoiViTi_hh
