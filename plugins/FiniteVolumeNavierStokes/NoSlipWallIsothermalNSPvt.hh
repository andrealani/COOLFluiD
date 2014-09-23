#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSPvt_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

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
class NoSlipWallIsothermalNSPvt : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallIsothermalNSPvt(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NoSlipWallIsothermalNSPvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

 protected:
  
  /// reference temperature
  CFreal m_refTemp;
  
  /// temperature at the wall
  CFreal m_wallTemp;
  
}; // end of class NoSlipWallIsothermalNSPvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSPvt_hh
