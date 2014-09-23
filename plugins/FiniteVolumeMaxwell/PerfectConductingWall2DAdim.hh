#ifndef COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2DAdim_hh
#define COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2DAdim_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell2DVarSetAdim;
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
class PerfectConductingWall2DAdim : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  PerfectConductingWall2DAdim(const std::string& name);
  
  /**
   * Default destructor
   */
  ~PerfectConductingWall2DAdim();
  
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
  Common::SafePtr<Physics::Maxwell::Maxwell2DVarSetAdim> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;
    
}; // end of class PerfectConductingWall2DAdim

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PerfectConductingWall2DAdim_hh
