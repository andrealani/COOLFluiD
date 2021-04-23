#ifndef COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaPuvt_hh
#define COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/NavierStokesKLogOmegaVarSet.hh"
#include "KOmega/NavierStokesKLogOmegaPvt.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables and Gamma-Alpha trasition model coupled with K-LogOmega turbulence model
   *
   * @author Ray Vandenhoeck
   */

class NavierStokes3DGammaAlphaPuvt : public 
KOmega::NavierStokesKLogOmegaPvt<KOmega::NavierStokesKLogOmegaVarSet<NavierStokes::NavierStokes3DVarSet, 0> > {
public: // classes
  
  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DGammaAlphaPuvt(const std::string& name,
			  Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DGammaAlphaPuvt();
  
  virtual void setComposition(const RealVector& state,
			      const bool isPerturb,
			      const CFuint iVar)
  {
    _isPerturb = isPerturb;
    _iPerturbVar = iVar;
  }
  
  /**
   * Get number of turbulent variables
   */
  virtual CFuint getNbTurbVars() const 
  {
    const CFuint nbTurbVars = _eulerModel->getNbScalarVars(0);
    cf_assert(nbTurbVars == 4);
    return nbTurbVars;
  }
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);
  
  /// Compute the stress tensor
  virtual void computeStressTensor(const RealVector& state,
				   const std::vector<RealVector*>& gradients,
				   const CFreal& radius); 

private :
  
  CFreal _unperturbedFluxGa;
  CFreal _unperturbedFluxAlpha;

}; // end of class NavierStokes3DGammaAlphaPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaPuvt_hh
