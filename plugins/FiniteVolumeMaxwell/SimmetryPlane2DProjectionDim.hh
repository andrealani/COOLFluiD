#ifndef COOLFluiD_Numerics_FiniteVolume_SimmetryPlane2DProjectionDim_hh
#define COOLFluiD_Numerics_FiniteVolume_SimmetryPlane2DProjectionDim_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell2DProjectionAdimVarSet;
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
class SimmetryPlane2DProjectionDim : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  SimmetryPlane2DProjectionDim(const std::string& name);
  
  /**
   * Default destructor
   */
  ~SimmetryPlane2DProjectionDim();
   
//   /**
//    * Defines the Config Option's of this class
//    * @param options a OptionList where to add the Option's
//    */
//   static void defineConfigOptions(Config::OptionList& options); 
//   
//   /**
//    * Configure the object
//    */  
//   virtual void configure ( Config::ConfigArgs& args );
  
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
  Common::SafePtr<Physics::Maxwell::Maxwell2DProjectionAdimVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// Electrical conductivity
  CFreal _electricalConductivity;
  
  /// Electrical conductivity
  CFreal _Ez0;
    
}; // end of class SimmetryPlane2DProjectionDim

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SimmetryPlane2DProjectionDim_hh
