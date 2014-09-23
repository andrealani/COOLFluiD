#ifndef COOLFluiD_Numerics_FluctSplit_WeakSubInletEuler2DConsRingleb_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSubInletEuler2DConsRingleb_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC2D.hh"
#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak sub inlet bc for Euler2D
 *
 * @author Andrea Lani
 *
 *
 *
 */
class WeakSubInletEuler2DConsRingleb : public WeakBC2D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSubInletEuler2DConsRingleb(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSubInletEuler2DConsRingleb();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

 protected:

/**
   * Set the state vector in the ghost State's
   */
  void setGhostState(const Framework::State& state,
		     Framework::State& gstate);

  virtual CFreal getAngle(const Framework::State& state);
 private:

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// flow angle
  CFreal _angle;
  /// total temperature
  CFreal m_tTotal;

  /// total pressure
  CFreal m_mach;

CFreal _pTotal;
}; // end of class WeakSubInletEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSubInletEuler2DCons_hh
