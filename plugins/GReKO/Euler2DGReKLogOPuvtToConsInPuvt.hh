#ifndef COOLFluiD_Physics_GReKO_Euler2DGReKLogOPuvtToConsInPuvt_hh
#define COOLFluiD_Physics_GReKO_Euler2DGReKLogOPuvtToConsInPuvt_hh

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
 * This class represents a transformer of variables from primitive [p u v T K Omega Gamma Retheta]
 * to conservative [rho rhoU rhoV rhoE rhoK rhoOmega rhoGamma rhoRetheta]  starting from primitive variables
 *
 * @author Khalil Bensassi 
 * @author Ray Vandenhoeck
 *
 */
class Euler2DGReKLogOPuvtToConsInPuvt : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGReKLogOTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DGReKLogOPuvtToConsInPuvt
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DGReKLogOPuvtToConsInPuvt();

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

}; // end of class Euler2DGReKLogOPuvtToConsInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_Euler2DGReKLogOPuvtToConsInPuvt_hh
