#ifndef COOLFluiD_Numerics_FiniteVolume_ICPInductionBC_hh
#define COOLFluiD_Numerics_FiniteVolume_ICPInductionBC_hh

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
class ICPInductionBC : public BASE {

public: 
  
  /**
   * Constructor
   */
  ICPInductionBC(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ICPInductionBC();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

}; // end of class ICPInductionBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICPInductionBC.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ICPInductionBC_hh
