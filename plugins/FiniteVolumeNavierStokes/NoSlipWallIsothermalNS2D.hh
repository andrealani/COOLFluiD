#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNS2D_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNS2D_hh

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
class NoSlipWallIsothermalNS2D : public FVMCC_BC { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  NoSlipWallIsothermalNS2D(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNS2D();
  
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
  
  /// temperature at the wall
  CFreal _wallTemp;
  
  /// temporary node
  RealVector _tempNode;
  
  /// temporary middle node
  RealVector _midNode;
  
  /// temporary ghost node
  RealVector _tempGhostNode;
  
}; // end of class NoSlipWallIsothermalNS2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNS2D_hh
