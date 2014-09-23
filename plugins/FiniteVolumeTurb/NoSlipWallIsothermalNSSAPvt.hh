#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSSAPvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSSAPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSPvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the mirror bc
   *
   * @author Andrea Lani
   *
   */
class NoSlipWallIsothermalNSSAPvt : public NoSlipWallIsothermalNSPvt {

public:

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSSAPvt(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallIsothermalNSSAPvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Add the value on the boundary for K in the turbulence model
   */
  void setGhostState(Framework::GeometricEntity *const face);


}; // end of class NoSlipWallIsothermalNSSAPvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSSAPvt_hh
