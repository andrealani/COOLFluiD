#ifndef COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall3DProjectionDim_hh
#define COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall3DProjectionDim_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell3DProjectionVarSet;
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
class PerfectConductingWall3DProjectionDim : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  PerfectConductingWall3DProjectionDim(const std::string& name);
  
  /**
   * Default destructor
   */
  ~PerfectConductingWall3DProjectionDim();
  
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
  Common::SafePtr<Physics::Maxwell::Maxwell3DProjectionVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;
    
}; // end of class PerfectConductingWall3DProjectionDim

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall3DProjectionDim_hh
