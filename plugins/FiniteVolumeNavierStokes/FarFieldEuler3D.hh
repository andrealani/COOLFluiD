#ifndef COOLFluiD_Numerics_FiniteVolume_FarFieldEuler3D_hh
#define COOLFluiD_Numerics_FiniteVolume_FarFieldEuler3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/FarFieldEuler2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to implement farfield boundary condition
 *
 * @author Andrea Lani
 *
 */
class FarFieldEuler3D : public FarFieldEuler2D {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FarFieldEuler3D(const std::string& name);

  /**
   * Default destructor
   */
  ~FarFieldEuler3D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

protected: // data

  /// z velocity
  CFreal                                   _wInf;

}; // end of class FarFieldEuler3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FarFieldEuler3D_hh
