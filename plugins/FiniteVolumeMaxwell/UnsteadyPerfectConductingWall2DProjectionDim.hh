#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadyPerfectConductingWall2DProjectionDim_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadyPerfectConductingWall2DProjectionDim_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Perfect Conducting wall boundary condition with projection scheme command 
   * 
   * @author Alejandro Alvarez Laguna
   *
   */
class UnsteadyPerfectConductingWall2DProjectionDim : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  UnsteadyPerfectConductingWall2DProjectionDim(const std::string& name);
  
  /**
   * Default destructor
   */
  ~UnsteadyPerfectConductingWall2DProjectionDim();
  
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
  Common::SafePtr<Physics::Maxwell::Maxwell2DProjectionVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;
    
}; // end of class UnsteadyPerfectConductingWall2DProjectionDim

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadyPerfectConductingWall2DProjectionDim_hh
