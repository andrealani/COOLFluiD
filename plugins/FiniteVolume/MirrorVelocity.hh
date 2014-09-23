#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorVelocity_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorVelocity_hh

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

class MirrorVelocity  : public FVMCC_BC { 

public: 
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * PuvtLTEtructor
   */
  MirrorVelocity(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~MirrorVelocity();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

protected:
  
  /// array of flags telling if a variable is a velocity component
  std::valarray<bool> m_isVelocityComp;
  
}; // end of class MirrorVelocity

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorVelocity_hh
