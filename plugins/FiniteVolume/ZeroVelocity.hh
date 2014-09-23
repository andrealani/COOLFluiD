#ifndef COOLFluiD_Numerics_FiniteVolume_ZeroVelocity_hh
#define COOLFluiD_Numerics_FiniteVolume_ZeroVelocity_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command applying a mirror boundary condition
   * for the velocity components
   * 
   * @author Andrea Lani
   */

//////////////////////////////////////////////////////////////////////////////

class ZeroVelocity  : public FVMCC_BC { 

public: 
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * PuvtLTEtructor
   */
  ZeroVelocity(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ZeroVelocity();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

private:
  
  /// array of flags telling if a variable is a velocity component
  std::valarray<bool> _isVelocityComp;
  
  /// array storing the IDs of the variables corresponding to velocity components
  std::vector<CFuint> _velocityIDs;

}; // end of class ZeroVelocity

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ZeroVelocity_hh
