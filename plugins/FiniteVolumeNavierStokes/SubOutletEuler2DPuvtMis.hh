#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtMis_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtMis_hh

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
 * @author Radek Honzatko
 *
 */
class SubOutletEuler2DPuvtMis : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletEuler2DPuvtMis(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubOutletEuler2DPuvtMis();

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

  /// isoentropic Mach number
  CFreal M_is;
 
  /// total pressure at inlet
  CFreal Pt_inlet;

}; // end of class SubOutletEuler2DPuvtMis

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEuler2DPuvtMis_hh
