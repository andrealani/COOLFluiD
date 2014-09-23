#ifndef COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler3DConsHO_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler3DConsHO_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC3DHO.hh"
#include "NavierStokes/Euler3DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak sub inlet bc for Euler3D
 *
 * @author Andrea Lani
 *
 *
 *
 */
class WeakSubOutletEuler3DConsHO : public WeakBC3DHO {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSubOutletEuler3DConsHO(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSubOutletEuler3DConsHO();

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
  Common::SelfRegistPtr<Physics::NavierStokes::Euler3DCons> _varSet;

  /// total temperature
  CFreal _pressure;

}; // end of class WeakSubOutletEuler3DConsHO

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler3DCons_hh
