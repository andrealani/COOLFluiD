#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletInterpYiVTTv_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletInterpYiVTTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/SuperInletInterp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////
        
 /**
  * This class implements a subsonic inlet that interpolates data from a file
  * for providing mass fractions, velocity, temperatures. Pressure is 
  * extrapolated from inside.
  *
  * @author Andrea Lani
  *
  */   
class SubInletInterpYiVTTv : public SuperInletInterp {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SubInletInterpYiVTTv(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubInletInterpYiVTTv();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
protected: // data
  
  /// pointer to the physical chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// species molar masses
  RealVector m_mmasses;
   
  /// mass fractions, roto-translational and vibrational temperatures in the ghost state
  RealVector m_ghostYiTTv;
  
  /// mass fractions, roto-translational and vibrational temperatures in the inner state
  RealVector m_innerYiTTv;
  
  /// mass fractions, roto-translational and vibrational temperatures on the boundary
  RealVector m_boundYiTTv;
  
  /// minimum values to fix for Yi, T, Tv
  RealVector m_minYiTTv;
  
  /// inlet values
  std::vector<CFreal> m_yiVTTv;
  
  /// blowing velocity
  CFreal m_blowVelocity;

  /// Vatsalya: Additional variables like phi
  CFreal m_addVar;
  
}; // end of class SubInletInterpYiVTTv

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletInterpYiVTTv_hh
