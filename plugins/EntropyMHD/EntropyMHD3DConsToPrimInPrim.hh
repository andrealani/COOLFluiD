#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHD3DConsToPrimInPrim_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHD3DConsToPrimInPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Jacobian matrix dW/dU (Cons->Prim) evaluated in Prim variables.
 *
 * EntropyMHD Cons = (rho, rho*u, rho*v, rho*w, Bx, By, Bz, sigma, phi)
 * EntropyMHD Prim = (rho,  u,  v,  w, Bx, By, Bz, sigma, phi)
 *
 * Only rows 1-3 differ from identity (velocity <-> momentum).
 * Row 7 (sigma) is identity; sigma is unchanged between Cons and Prim.
 *
 * @author Rayan Dhib
 */
class EntropyMHD3DConsToPrimInPrim :
        public Framework::VarSetMatrixTransformer {
public:

  EntropyMHD3DConsToPrimInPrim(Common::SafePtr<Framework::PhysicalModelImpl> model);

  ~EntropyMHD3DConsToPrimInPrim();

  /**
   * Set the transformation matrix from a given state (in Prim variables)
   */
  void setMatrix(const RealVector& state);

private:

  bool getIsIdentityTransformation() const
  {
    return false;
  }

}; // end of class EntropyMHD3DConsToPrimInPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHD3DConsToPrimInPrim_hh
