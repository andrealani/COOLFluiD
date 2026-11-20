#ifndef COOLFluiD_Physics_MHD_MHDProjectionDiffPrimRes_hh
#define COOLFluiD_Physics_MHD_MHDProjectionDiffPrimRes_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the diffusive terms of resistive and viscous MHD in
 * primitive variables (rho, u, v, w, Bx, By, Bz, p, phi), with a uniform
 * resistivity eta and a uniform kinematic viscosity nu, in 2D and 3D.
 *
 * The diffusive fluxes through a face with normal n are:
 *  - momentum: F(rho v_i) = nu rho grad(v_i) . n
 *  - induction: F(B_i) = eta sum_j (d_j B_i - d_i B_j) n_j, i.e. -curl(eta curl B)
 *  - total energy: F(E) = sum_i B_i F(B_i) + sum_i v_i F(rho v_i)
 * where d_j is the derivative along direction j (zero along z in 2D).
 *
 * @author Rayan Dhib
 */
class MHDProjectionDiffPrimRes : public MHDProjectionDiffVarSet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  MHDProjectionDiffPrimRes(const std::string& name,
                           Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~MHDProjectionDiffPrimRes();

  /**
   * Set up private data (reads the space dimension)
   */
  virtual void setup();

  /**
   * Set the composition
   */
  virtual void setComposition(const RealVector& state,
                              const bool isPerturb,
                              const CFuint iVar);

  /**
   * Set the quantities needed to compute gradients starting from the states
   */
  virtual void setGradientVars(const std::vector<RealVector*>& states,
                               RealMatrix& values,
                               const CFuint stateSize);

  /**
   * Compute required gradients starting from the gradients of the states
   */
  virtual void setGradientVarGradients(const std::vector<RealVector*>& states,
                                       const std::vector< std::vector<RealVector*> >& stateGradients,
                                       std::vector< std::vector<RealVector*> >& gradVarGradients,
                                       const CFuint stateSize);

  /**
   * Compute the gradients of the states starting from gradient variable gradients
   */
  virtual void setStateGradients(const std::vector<RealVector*>& states,
                                 const std::vector< std::vector<RealVector*> >& gradVarGradients,
                                 std::vector< std::vector<RealVector*> >& stateGradients,
                                 const CFuint stateSize);

  /// Get the diffusive flux Jacobian
  virtual void computeFluxJacobian(const RealVector& state,
                                   const RealVector& gradientJacob,
                                   const RealVector& normal,
                                   const CFreal& radius,
                                   RealMatrix& fluxJacob);

  /**
   * Get the adimensional dynamic viscosity
   */
  virtual CFreal getDynViscosity(const RealVector& state,
                                 const std::vector<RealVector*>& gradients);

  /**
   * Get the adimensional density
   */
  virtual CFreal getDensity(const RealVector& state);

  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
                                      const CFreal& dynViscosity)
  {
    return 0.;
  }

  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

  /**
   * Get the diffusive flux vector (not implemented)
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);

  /**
   * Get the current update diffusive coefficient for CFL
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS)
  {
    return std::max(_eta, _nu);
  }

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

private:

  /// magnetic resistivity
  CFreal _eta;

  /// kinematic viscosity
  CFreal _nu;

  /// number of space dimensions (2 or 3)
  CFuint _dim;

}; // end of class MHDProjectionDiffPrimRes

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionDiffPrimRes_hh
