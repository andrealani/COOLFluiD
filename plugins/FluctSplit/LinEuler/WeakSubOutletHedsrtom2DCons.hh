#ifndef COOLFluiD_Numerics_FluctSplit_WeakSubOutletHedstrom2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSubOutletHedstrom2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakBC2D.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "LinEuler/LinEuler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak sub outlet bc for LinEuler2D
 *
 * @author Lilla Koloszar
 *
 *
 *
 */
class WeakSubOutletHedstrom2DCons : public WeakBC2D {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSubOutletHedstrom2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSubOutletHedstrom2DCons();

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
  Common::SelfRegistPtr<Physics::LinearizedEuler::LinEuler2DCons> _varSet;

  /// normals of the nodes
  std::vector< std::vector<RealVector> > _bcNormals;

}; // end of class WeakSubOutletHedstrom2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSubOutletHedstrom2DCons_hh
