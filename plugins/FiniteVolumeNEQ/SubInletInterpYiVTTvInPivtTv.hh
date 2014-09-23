#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletInterpYiVTTvInPivtTv_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletInterpYiVTTvInPivtTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNEQ/SubInletInterpYiVTTv.hh"

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
class SubInletInterpYiVTTvInPivtTv : public SubInletInterpYiVTTv {
public:
  
  /**
   * Constructor
   */
  SubInletInterpYiVTTvInPivtTv(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubInletInterpYiVTTvInPivtTv();

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
  
  /// array of partial densities
  RealVector m_rhoi;
  
  /// array of partial densities on the boundary
  RealVector m_rhoiB;
  
}; // end of class SubInletInterpYiVTTvInPivtTv

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletInterpYiVTTvInPivtTv_hh
