#ifndef COOLFluiD_Numerics_FiniteVolume_ICPMirror_hh
#define COOLFluiD_Numerics_FiniteVolume_ICPMirror_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolumeICP {
          
//////////////////////////////////////////////////////////////////////////////
      
      /**
       * This class represents a command that applies the BC for the 
       * induction equation in ICP simulation
       * 
       * @author Andrea Lani
       * @author Radek Honzatko
       */
template <class BASE>
class ICPMirror : public BASE {

public: 
  
  /**
   * Constructor
   */
  ICPMirror(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ICPMirror();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

}; // end of class ICPMirror

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICPMirror.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ICPMirror_hh
