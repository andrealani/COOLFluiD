#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWall2DProjectionDim_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWall2DProjectionDim_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Perfect Conducting wall boundary condition with projection scheme command 
   * 
   * @author Alejandro Alvarez Laguna
   *
   */
class NoSlipWall2DProjectionDim : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  NoSlipWall2DProjectionDim(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NoSlipWall2DProjectionDim();
   
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options); 

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
  Common::SafePtr<Physics::Maxwell::Maxwell2DProjectionVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// Electrical conductivity
  CFreal _electricalConductivity;
  
  /// Electrical conductivity
  CFreal _Ez0;
  
  ///introduced non induced electromagnetic field
  std::vector<CFreal> _nonInducedEMField;
    
}; // end of class NoSlipWall2DProjectionDim

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWall2DProjectionDim_hh
