#ifndef COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOSSTPuvt_hh
#define COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOSSTPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes3DGReKLogOBSLPuvt.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables and Menter's K-LogOmega SST turbulence model
   *
   * @author Khalil Bensassi 
   * @author Ray Vandenhoeck
   */

class NavierStokes3DGReKLogOSSTPuvt : public NavierStokes3DGReKLogOBSLPuvt {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGReKLogOTerm;

  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DGReKLogOSSTPuvt(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DGReKLogOSSTPuvt();

  /**
   * Set the coefficients for the model
   */
  virtual void setModelCoefficients();

  /**
   * Get the adimensional turbulent dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getTurbDynViscosityFromGradientVars(const RealVector& state, const std::vector<RealVector*>& gradients);

}; // end of class NavierStokes3DGReKLogOSSTPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOSSTPuvt_hh
