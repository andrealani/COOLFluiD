#ifndef COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2DHalf_hh
#define COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2DHalf_hh

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
class PerfectConductingWall2DHalf : public FVMCC_BC {

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  PerfectConductingWall2DHalf(const std::string& name);
  
  /**
   * Default destructor
   */
  ~PerfectConductingWall2DHalf();
  
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
  
    
}; // end of class PerfectConductingWall2DHalf

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2DHalf_hh
