#ifndef COOLFluiD_Physics_KOmega_Euler3DKLogOmegaConsToPuvtInPuvt_hh
#define COOLFluiD_Physics_KOmega_Euler3DKLogOmegaConsToPuvtInPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative [rho rhoU rhoV rhoW rhoE rhoK rhoOmega]
 * to primitive [p u v w T K Omega] starting from primitive variables
 *
 * @author Thomas Wuilbaut
 * @author Milan Zaloudek
 * @author Ray Vandenhoeck
 *
 */
class Euler3DKLogOmegaConsToPuvtInPuvt : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerKLogOmegaTerm;

  /**
   * Default constructor without arguments
   */
  Euler3DKLogOmegaConsToPuvtInPuvt
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DKLogOmegaConsToPuvtInPuvt();

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
  Common::SafePtr<EulerKLogOmegaTerm> _model;

}; // end of class Euler3DKLogOmegaConsToPuvtInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_Euler3DKLogOmegaConsToPuvtInPuvt_hh
