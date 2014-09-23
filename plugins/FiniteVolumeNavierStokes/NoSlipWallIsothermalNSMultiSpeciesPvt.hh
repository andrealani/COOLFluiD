#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSMultiSpeciesPvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSMultiSpeciesPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSPvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the Isothermal wall
   * bc with LTE Demixing
   * 
   * @author Janos Molnar 
   * @author Andrea Lani
   *
   */
template <class PMODEL>
class NoSlipWallIsothermalNSMultiSpeciesPvt : public NoSlipWallIsothermalNSPvt {

public: 
  
  /**
   * Constructor
   */
  NoSlipWallIsothermalNSMultiSpeciesPvt(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNSMultiSpeciesPvt();
  
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
  
  /// number of chemical variables (elements, species, ...)
  CFuint _nbChemVars;
  
}; // end of class NoSlipWallIsothermalNSMultiSpeciesPvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSMultiSpeciesPvt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSMultiSpeciesPvt_hh
