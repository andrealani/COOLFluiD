#ifndef COOLFluiD_Numerics_FluctSplit_WeakSubInletEuler3DCons_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSubInletEuler3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC3D.hh"
#include "NavierStokes/Euler3DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak sub inlet bc for Euler3D
 *
 * @author Andrea Lani
 * @author Fabio Pinna
 */
class WeakSubInletEuler3DCons : public WeakBC3D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSubInletEuler3DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSubInletEuler3DCons();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

 protected:

  /**
   * Set the additional flux and the jacobian of the fluxes
   */
  void setGhostState(const Framework::State& state,
		     Framework::State& gstate);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

private:

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler3DCons> m_varSet;

  /// temporary storage for the jacobian of the gstate
  RealMatrix m_dUinDu;

  /// storage for the identity matrix (only the diagonal is stored)
  RealVector m_identity;

  /// total temperature
  CFreal m_tTotal;

  /// total pressure
  CFreal m_pTotal;

  /// flow angle (Azimuth)
  CFreal m_angleAzi;
  
  /// flow angle (Zenith)
  CFreal m_angleZen;

}; // end of class WeakSubInletEuler3DConsImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSubInletEuler3DConsImpl_hh
