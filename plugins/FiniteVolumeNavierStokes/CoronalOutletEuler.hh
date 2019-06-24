#ifndef COOLFluiD_Numerics_FiniteVolume_CoronalOutletEuler_hh
#define COOLFluiD_Numerics_FiniteVolume_CoronalOutletEuler_hh

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
class CoronalOutletEuler : public FVMCC_BC {
public: 
  
  /**
   * Constructor
   */
  CoronalOutletEuler(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~CoronalOutletEuler();
  
  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
}; // end of class CoronalOutletEuler

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CoronalOutletEuler_hh
