#ifndef COOLFluiD_Numerics_FiniteVolume_NullBC_hh
#define COOLFluiD_Numerics_FiniteVolume_NullBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class applies a Null boundary condition
   * 
   * @author Andrea Lani
   */

//////////////////////////////////////////////////////////////////////////////

class NullBC  : public FVMCC_BC { 

public: 

  /**
   * PuvtLTEtructor
   */
  NullBC(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NullBC();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

}; // end of class NullBC

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NullBC_hh
