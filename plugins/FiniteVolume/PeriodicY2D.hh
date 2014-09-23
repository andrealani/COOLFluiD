#ifndef COOLFluiD_Numerics_FiniteVolume_PeriodicY2D_hh
#define COOLFluiD_Numerics_FiniteVolume_PeriodicY2D_hh

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
   * This class represents a periodic boundary condition command for two 
   * topological regions in 2D parallel to y-axis for CellCenterFVM schemes
   * for serial simulations
   *
   * @author Mehmet Sarp Yalim
   *
   */
class PeriodicY2D : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  PeriodicY2D(const std::string& name);

  /**
   * Default destructor
   */
  ~PeriodicY2D();
  
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
  
  /// vector of boundary states corresponding to the PeriodicY faces
  std::vector<Framework::State*> _boundaryStates;

  /// map that maps global topological region face ID to a local one in the topological region set 
  Common::CFMap<CFuint,CFuint> _globalToLocalTRSFaceID; 
  
  /// threshold
  CFreal _threshold;
  
}; // end of class PeriodicY2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PeriodicY2D_hh
