#ifndef COOLFluiD_Physics_SA_Euler2DSAPuvtToConsInPuvt_hh
#define COOLFluiD_Physics_SA_Euler2DSAPuvtToConsInPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive [p u v T rhoK]
 * to conservative [rho rhoU rhoV rhoE rhoK] starting from primitive variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DSAPuvtToConsInPuvt : public Framework::VarSetMatrixTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerSATerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DSAPuvtToConsInPuvt
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DSAPuvtToConsInPuvt();

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
  Common::SafePtr<EulerSATerm> _model;

}; // end of class Euler2DSAPuvtToConsInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_Euler2DSAPuvtToConsInPuvt_hh
