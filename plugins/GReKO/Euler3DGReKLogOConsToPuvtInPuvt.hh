#ifndef COOLFluiD_Physics_GReKO_Euler3DGReKLogOConsToPuvtInPuvt_hh
#define COOLFluiD_Physics_GReKO_Euler3DGReKLogOConsToPuvtInPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative [rho rhoU rhoV rhoE rhoK rhoOmega rhoGamma rhoReTheta]
 * to primitive [p u v w T K Omega Gamma ReTheta] starting from primitive variables
 *
 * @author Khalil Bensassi
 * @author Ray Vandenhoeck
 *
 */
class Euler3DGReKLogOConsToPuvtInPuvt : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGReKLogOTerm;

  /**
   * Default constructor without arguments
   */
  Euler3DGReKLogOConsToPuvtInPuvt
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DGReKLogOConsToPuvtInPuvt();

  /**
   * Set the transformation matrix from a given state
   */
  void setMatrix(const RealVector& state);

private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerGReKLogOTerm> _model;

}; // end of class Euler3DGReKLogOConsToPuvtInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_Euler3DGReKLogOConsToPuvtInPuvt_hh
