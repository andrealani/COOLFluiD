#ifndef COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar3D_hh
#define COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/FarFieldEuler3D.hh"

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
class FarFieldEulerChar3D : public FarFieldEuler3D {
public:

  /**
   * Constructor
   */
  FarFieldEulerChar3D(const std::string& name);

  /**
   * Default destructor
   */
  ~FarFieldEulerChar3D();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

}; // end of class FarFieldEulerChar3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FarFieldEulerChar3D_hh
