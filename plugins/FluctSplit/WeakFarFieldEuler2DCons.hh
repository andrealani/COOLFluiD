#ifndef COOLFluiD_Numerics_FluctSplit_WeakFarFieldEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_WeakFarFieldEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC2D.hh"
#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak bc for Euler 2D
 * It checks if the flow is comming in or out
 * And apply the correponding condition
 *
 * @author Nadege Villedieu
 *
 *
 *
 */
class WeakFarFieldEuler2DCons : public WeakBC2D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakFarFieldEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakFarFieldEuler2DCons();

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

  
  /// total temperature
  CFreal _tTotal;

  /// total pressure
  CFreal _pTotal;

  /// flow angle
  CFreal _angle;

}; // end of class WeakFarFieldEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakFarFieldEuler2DCons_hh
