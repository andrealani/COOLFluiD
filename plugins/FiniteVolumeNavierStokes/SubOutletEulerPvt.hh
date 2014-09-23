#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletEulerPvt_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletEulerPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic outlet command
 *
 * @author Andrea Lani
 *
 */
template <class VARSET>
class SubOutletEulerPvt : public FVMCC_BC {
public:

  typedef typename VARSET::EULERSET EULERSET;

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletEulerPvt(const std::string& name);

  /**
   * Default destructor
   */
  ~SubOutletEulerPvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private: // data

  /// static pressure or pressure variation
  CFreal m_pressure;

}; // end of class SubOutletEulerPvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubOutletEulerPvt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEulerPvt_hh
