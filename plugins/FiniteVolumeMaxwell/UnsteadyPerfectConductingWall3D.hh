#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadyPerfectConductingWall3D_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadyPerfectConductingWall3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell3DVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Perfect Conducting wall boundary condition command 
   * 
   * @author Alejandro Alvarez Laguna
   *
   */
class UnsteadyPerfectConductingWall3D : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  UnsteadyPerfectConductingWall3D(const std::string& name);
  
  /**
   * Default destructor
   */
  ~UnsteadyPerfectConductingWall3D();
  
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
  Common::SafePtr<Physics::Maxwell::Maxwell3DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;
    
}; // end of class UnsteadyPerfectConductingWall3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadyPerfectConductingWall3D_hh
