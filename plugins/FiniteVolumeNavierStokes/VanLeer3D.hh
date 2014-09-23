#ifndef COOLFluiD_Numerics_FiniteVolume_VanLeer3D_hh
#define COOLFluiD_Numerics_FiniteVolume_VanLeer3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/VanLeer2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Van Leer flux splitting corresponding
 * to the Euler physical model 3D (in conservative variables)
 *
 * @author Janos Molnar
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class VanLeer3D : public VanLeer2D<UPDATEVAR> {
public: // classes

  /**
   * Constructor
   */
  VanLeer3D(const std::string& name);

  /**
   * Default destructor
   */
  ~VanLeer3D();
  
  /**
   * Compute the flux in the current face
   */
  virtual void compute(RealVector& result);
  
}; // end of class VanLeer3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "VanLeer3D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_VanLeer3D_hh
