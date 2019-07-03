#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInletPhotosphereHydrodynamic_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInletPhotosphereHydrodynamic_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/SuperInletProjection.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Supersonic inlet imposing photospheric 
   * conditions
   *
   * @author Andrea Lani
   * @author Peter Leitner
   *
   */
class SuperInletPhotosphereHydrodynamic : public SuperInletProjection {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperInletPhotosphereHydrodynamic(const std::string& name);

  /**
   * Default destructor
   */
  ~SuperInletPhotosphereHydrodynamic();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);
  
 private: // data
  
}; // end of class SuperInletPhotosphereHydrodynamic

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInletPhotosphereHydrodynamic_hh
