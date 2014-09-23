#ifndef COOLFluiD_Numerics_FiniteVolume_ShiftedPeriodicX2D_hh
#define COOLFluiD_Numerics_FiniteVolume_ShiftedPeriodicX2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a shifted periodic boundary condition command for two 
   * topological regions in 2D parallel to x-axis for CellCenterFVM schemes
   * for serial simulations
   *
   * @author Andrea Lani
   * @author Mehmet Sarp Yalim
   *
   */
class ShiftedPeriodicX2D : public FVMCC_BC {
public:

  /**
   * Constructor
   */
  ShiftedPeriodicX2D(const std::string& name);

  /**
   * Default destructor
   */
  ~ShiftedPeriodicX2D();
  
  /**
   * Set up private data
   */
  void setup ();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private: //data

  /// vector of connectivity of boundary faces corresponding to periodicity of topological regions
  std::vector<CFuint> _periodicFaceID;
  
  /// vector of boundary states corresponding to the ShiftedPeriodicX faces in sequence from left to right of the domain
  std::vector<Framework::State*> _boundaryStatesInSequence;

  /// map that maps global topological region face ID to a local one in the topological region set in sequence 
  Common::CFMap<CFuint,CFuint> _localTRSFaceIDInSequence; 

}; // end of class ShiftedPeriodicX2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ShiftedPeriodicX2D_hh
