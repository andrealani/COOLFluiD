#ifndef COOLFluiD_Physics_SA_NavierStokesSAVarSet_hh
#define COOLFluiD_Physics_SA_NavierStokesSAVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSet.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model for primitive
   * variables and SA turbulence model
   *
   * @author Joao Pinto
   * @author Thomas Wuilbaut
   * @author Khalil Bensassi
   */
template <typename BASE, int SGROUP>
class NavierStokesSAVarSet : public NavierStokes::NavierStokesTurbVarSet<BASE, SGROUP> {
  
public: // classes
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerSATerm;
  
  /**
   * Constructor
   */
  NavierStokesSAVarSet(const std::string& name,
		       Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesSAVarSet();
  
  /**
   * Get the adimensional dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients);
    
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
    
protected:

  /// convective model
  Common::SafePtr<EulerSATerm> _eulerModel;
  
  /// flag to know if we are computing the jacobian
  bool _isPerturb;

  /// flag to know if which variable we are perturbing
  CFuint _iPerturbVar;

  /// storage of the unperturbed flux for the first equation
  CFreal _unperturbedFluxNutil;
    
}; // end of class NavierStokesSAVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesSAVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_NavierStokesSAVarSet_hh
