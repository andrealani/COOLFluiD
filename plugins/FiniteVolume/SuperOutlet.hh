#ifndef COOLFluiD_Numerics_FiniteVolume_SuperOutlet_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperOutlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic inlet command 
   * 
   * @author Andrea Lani
   *
   */
class SuperOutlet : public FVMCC_BC {
public: 
  
  /**
   * Constructor
   */
  SuperOutlet(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~SuperOutlet();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
}; // end of class SuperOutlet

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperOutlet_hh
