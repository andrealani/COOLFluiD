#ifndef COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2D_hh
#define COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell2DVarSet;
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
class PerfectConductingWall2D : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  PerfectConductingWall2D(const std::string& name);
  
  /**
   * Default destructor
   */
  ~PerfectConductingWall2D();
  
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
  Common::SafePtr<Physics::Maxwell::Maxwell2DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;
    
}; // end of class PerfectConductingWall2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2D_hh
