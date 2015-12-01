#ifndef COOLFluiD_Physics_KOmega_NavierStokesKOmegaVarSet_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKOmegaVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSet.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and K-Omega turbulence model
   *
   * @author Thomas Wuilbaut
   * @author Khalil Bensassi
   */
template <typename BASE, int SGROUP>
class NavierStokesKOmegaVarSet : public NavierStokes::NavierStokesTurbVarSet<BASE, SGROUP> {

public: // classes
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerKOmegaTerm;
  
  /**
   * Constructor
   */
  NavierStokesKOmegaVarSet(const std::string& name,
			   Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesKOmegaVarSet();

  /**
   * Get the adimensional dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients);
  
  /**
   * Get the adimensional turbulent dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getTurbDynViscosityFromGradientVars(const RealVector& state, const std::vector<RealVector*>& gradients);
  
  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  virtual void setComposition(const RealVector& state,
			      const bool isPerturb,
			      const CFuint iVar)
  {
    _isPerturb = isPerturb;
    _iPerturbVar = iVar;
  }
  
  /**
   * set up private data
   */
  virtual void setup();
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);
  
  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);
  
  /**
   * Get the heat flux
   */
  virtual CFreal getHeatFlux(const RealVector& state,
                             const std::vector<RealVector*>& gradients,
                             const RealVector& normal);
  
  /**
   * Compute the blending coeficients needed for BSL and SST variants of the model
   */
  virtual void computeBlendingCoefFromGradientVars(const RealVector& state, const RealVector& gradK, const RealVector& gradOmega);
  
  /**
   * Return the blending coefficient F1 (or F2 in the SST case)
   */
  CFreal getBlendingCoefficientF1() {return _blendingCoef1;}
  
  /**
   * Return the coefficient sigmaOmega2
   */
  CFreal getSigmaOmega2() {return _sigmaOmega2;}
  
  /**
   * Get the coefficients for the model
   */
  CFreal getSigmaK()
  {
    return ((_blendingCoef1 * _sigmaK1) + ((1.-_blendingCoef1)*_sigmaK2));
  }
  
  CFreal getSigmaOmega()
  {
    return ((_blendingCoef1 * _sigmaOmega1) + ((1.-_blendingCoef1)*_sigmaOmega2));
  }
  
  CFreal getGammaCoef()
  {
    return ((_blendingCoef1 * _gamma1) + ((1.-_blendingCoef1)*_gamma2));
  }

  /**
   * Get the compressible value of beta
   * using the model of Wilcox (Turbulence Modeling for CFD, DCW Industries,Inc.)
   */
  CFreal getBeta(const RealVector& state)
  {
    const CFreal betaIncomp = getIncompBeta();
    const CFreal betaStarIncomp = getIncompBetaStar();
    
    const CFreal Mt0 = 0.25;
    const CFreal ksiStar = 1.5;
    const CFreal K = state[this->_kID];
    const CFreal A2 = getSqSoundSpeed(state);
    const CFreal Mt = 2.*K/A2;
    
    //Heaviside function
    CFreal H = 0.;
    if((Mt - Mt0) == 0.) H = 0.5;
    if((Mt - Mt0) > 0.) H = 1.;
    
    const CFreal FMt = (Mt*Mt - Mt0*Mt0) * H;
    return betaIncomp - (betaStarIncomp * ksiStar * FMt);
  }

  /**
   * Get the compressible value of betaStar
   * using the model of Wilcox (Turbulence Modeling for CFD, DCW Industries,Inc.)
   */
  CFreal getBetaStar(const RealVector& state)
  {
    const CFreal Mt0 = 0.25;
    const CFreal betaStarIncomp = getIncompBetaStar();
    const CFreal ksiStar = 1.5;
    const CFreal K = state[this->_kID];
    const CFreal A2 = getSqSoundSpeed(state);
    const CFreal Mt = 2.*K/A2;
    
    //Heaviside function
    CFreal H = 0.;
    if((Mt - Mt0) == 0.) H = 0.5;
    if((Mt - Mt0) > 0.) H = 1.;
    
    const CFreal FMt = (Mt*Mt - Mt0*Mt0) * H;
    return betaStarIncomp * (1. + (ksiStar * FMt));
  }
  
  /**
   * Get number of turbulent variables
   */
  virtual CFuint getNbTurbVars() const 
  {
    const CFuint nbTurbVars = _eulerModel->getNbScalarVars(SGROUP);
    cf_assert(nbTurbVars == 2);
    return nbTurbVars;
  }
  
protected :

  /**
   * Set the coefficients for the model
   */
  virtual void setModelCoefficients();
  
  /// Compute the square of the speed of sound
  virtual CFreal getSqSoundSpeed(const RealVector& state) = 0;
  
private:
  
  /**
   * Get the incompressible value of beta
   */
  CFreal getIncompBeta()
  {
    return ((_blendingCoef1 * _beta1) + ((1.-_blendingCoef1)*_beta2));
  }

  /**
   * Get the incompressible value of betaStar
   */
  CFreal getIncompBetaStar() {return _betaStar;}
  
protected :

  /// convective model
  Common::SafePtr<EulerKOmegaTerm> _eulerModel;
  
  //storage of the blending coefficient
  CFreal _blendingCoef1;

  //storage of the model coefficients
  CFreal _sigmaK1;
  CFreal _sigmaK2;
  CFreal _sigmaOmega1;
  CFreal _sigmaOmega2;
  CFreal _beta1;
  CFreal _beta2;
  CFreal _gamma1;
  CFreal _gamma2;
  CFreal _betaStar;
  CFreal _kappa;

  /// flag to know if we are computing the jacobian
  bool _isPerturb;

  /// flag to know if which variable we are perturbing
  CFuint _iPerturbVar;

  /// storage of the unperturbed flux for the first equation
  CFreal _unperturbedFluxK;

  /// storage of the unperturbed flux for the second equation
  CFreal _unperturbedFluxOmega;
  
}; // end of class NavierStokesKOmegaVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKOmegaVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKOmegaVarSet_hh
