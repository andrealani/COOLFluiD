#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerWeakSubInletEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerWeakSubInletEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "TwoLayerWeakBC2D.hh"
#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak sub inlet bc for Euler2D
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerWeakSubInletEuler2DCons : public TwoLayerWeakBC2D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  TwoLayerWeakSubInletEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~TwoLayerWeakSubInletEuler2DCons();

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
  CFreal _tTotal;

  /// total pressure
  CFreal _pTotal;

  /// flow angle
  CFreal _angle;

}; // end of class TwoLayerWeakSubInletEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerWeakSubInletEuler2DCons_hh
