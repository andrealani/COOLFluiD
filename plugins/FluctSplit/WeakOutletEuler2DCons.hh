#ifndef COOLFluiD_Numerics_FluctSplit_WeakOutletEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_WeakOutletEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC2D.hh"
#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak subsonic-supersonic outlet bc for Euler2D
 *
 * @author Andrea Lani
 */
class WeakOutletEuler2DCons : public WeakBC2D {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakOutletEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakOutletEuler2DCons();

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
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

 protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Set the state vector in the ghost State's
   */
  void setGhostState(const Framework::State& state,
		     Framework::State& gstate);

 private:

  /// socket for the isUpdated data
  Framework::DataSocketSink< bool> socket_isUpdated;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> m_varSet;

  /// total temperature
  CFreal m_totTemperature;

}; // end of class WeakOutletEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakOutletEuler2DCons_hh
