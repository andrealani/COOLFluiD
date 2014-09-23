#ifndef COOLFluiD_Physics_GReKO_NavierStokes2DGReKOSSTPuvt_hh
#define COOLFluiD_Physics_GReKO_NavierStokes2DGReKOSSTPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DGReKOBSLPuvt.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and Menter's K-Omega SST turbulence model
   *
   * @author Khalil Bensassi 
   */

class NavierStokes2DGReKOSSTPuvt : public NavierStokes2DGReKOBSLPuvt {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGReKOTerm;

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DGReKOSSTPuvt(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes2DGReKOSSTPuvt();

  /**
   * Set the coefficients for the model
   */
  virtual void setModelCoefficients();

  /**
   * Get the adimensional turbulent dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getTurbDynViscosityFromGradientVars(const RealVector& state, const std::vector<RealVector*>& gradients);

}; // end of class NavierStokes2DGReKOSSTPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokes2DGReKOSSTPuvt_hh
