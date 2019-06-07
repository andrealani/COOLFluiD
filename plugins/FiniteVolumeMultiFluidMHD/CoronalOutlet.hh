#ifndef COOLFluiD_Numerics_FiniteVolume_CoronalOutlet_hh
#define COOLFluiD_Numerics_FiniteVolume_CoronalOutlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic outlet command 
   * 
   * @author Andrea Lani
   * @author Peter Leitner
   *
   */
class CoronalOutlet : public FVMCC_BC {
public: 
  
  /**
   * Constructor
   */
  CoronalOutlet(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~CoronalOutlet();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
}; // end of class CoronalOutlet

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CoronalOutlet_hh
