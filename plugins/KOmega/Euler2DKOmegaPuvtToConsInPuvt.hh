#ifndef COOLFluiD_Physics_KOmega_Euler2DKOmegaPuvtToConsInPuvt_hh
#define COOLFluiD_Physics_KOmega_Euler2DKOmegaPuvtToConsInPuvt_hh

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
 * This class represents a transformer of variables from primitive [p u v T K Omega]
 * to conservative [rho rhoU rhoV rhoE rhoK rhoOmega] starting from primitive variables
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DKOmegaPuvtToConsInPuvt : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerKOmegaTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DKOmegaPuvtToConsInPuvt
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DKOmegaPuvtToConsInPuvt();

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
  Common::SafePtr<EulerKOmegaTerm> _model;

}; // end of class Euler2DKOmegaPuvtToConsInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_Euler2DKOmegaPuvtToConsInPuvt_hh
