#ifndef COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler2DCons_hh

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
class WeakSubOutletEuler2DCons : public WeakBC2D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSubOutletEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSubOutletEuler2DCons();

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

 private:

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// total temperature
  CFreal _pressure;

}; // end of class WeakSubOutletEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler2DCons_hh
