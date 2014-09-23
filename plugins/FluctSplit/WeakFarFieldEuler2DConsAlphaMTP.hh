#ifndef COOLFluiD_Numerics_FluctSplit_WeakFarFieldEuler2DConsAlphaMTP_hh
#define COOLFluiD_Numerics_FluctSplit_WeakFarFieldEuler2DConsAlphaMTP_hh

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
class WeakFarFieldEuler2DConsAlphaMTP : public WeakBC2D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakFarFieldEuler2DConsAlphaMTP(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakFarFieldEuler2DConsAlphaMTP();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get Flow Angle
   */
  CFreal getAlpha()
  {
    return _alpha;
  }


 protected:

  /**
   * Set the state vector in the ghost State's
   */
  void setGhostState(const Framework::State& state,
                     Framework::State& gstate);

 private:

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// flow angle in degrees
  CFreal _alpha;

  /// mach number
  CFreal _mach;

  /// temperature
  CFreal _temperature;

  /// pressure
  CFreal _pressure;

 }; // end of class WeakFarFieldEuler2DConsAlphaMTP

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakFarFieldEuler2DConsAlphaMTP_hh
