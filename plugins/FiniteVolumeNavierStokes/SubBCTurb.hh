#ifndef COOLFluiD_Numerics_FiniteVolume_SubBCTurb_hh
#define COOLFluiD_Numerics_FiniteVolume_SubBCTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a subsonic inlet condition for turbulent flows
  * the turbulent variables
  *
  * @author Andrea Lani
  */
template <typename BASE, typename UPDATEVS, bool INLET>
class SubBCTurb : public BASE {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubBCTurb(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubBCTurb();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

private: // data
  
  /// physical model var set
  Common::SafePtr<UPDATEVS> m_varSetTurb;
  
  /// physical model var set
  std::vector<CFreal>       m_turbVars;
  
}; // end of class SubBCTurb

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/SubBCTurb.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubBCTurb_hh
