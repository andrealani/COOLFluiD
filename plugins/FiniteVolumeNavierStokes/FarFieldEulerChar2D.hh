#ifndef COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar2D_hh
#define COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/FarFieldEuler2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to implement farfield boundary condition
 * using th characteristics
 *
 * @author Thomas Wuilbaut
 *
 */
class FarFieldEulerChar2D : public FarFieldEuler2D {
public:

  /**
   * Constructor
   */
  FarFieldEulerChar2D(const std::string& name);

  /**
   * Default destructor
   */
  ~FarFieldEulerChar2D();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

}; // end of class FarFieldEulerChar2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar2D_hh
