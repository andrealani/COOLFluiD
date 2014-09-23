#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorEuler2D_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the mirror bc 
   * 
   * @author Mehmet Sarp Yalim
   * @author Andrea Lani 
   *
   */
class MirrorEuler2D : public FVMCC_BC { 

public: 
  
  /**
   * Constructor
   */
  MirrorEuler2D(const std::string& name);
  
  /**
   * Default destructor
   */
  ~MirrorEuler2D();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private:
  
  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;
  
  /// physical model data
  RealVector _dataInnerState;
  
  /// physical model data
  RealVector _dataGhostState;
  
}; // end of class MirrorEuler2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorEuler2D_hh
