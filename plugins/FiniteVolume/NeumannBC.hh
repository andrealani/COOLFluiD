#ifndef COOLFluiD_Numerics_FiniteVolume_NeumannBC_hh
#define COOLFluiD_Numerics_FiniteVolume_NeumannBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Neumann boundary condition
   *
   * @author Andrea Lani
   *
   */
class NeumannBC : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NeumannBC(const std::string& name);

  /**
   * Default destructor
   */
  ~NeumannBC();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

}; // end of class NeumannBC

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NeumannBC_hh
