#ifndef COOLFluiD_Physics_GReKO_NavierStokes2DGReKOPuvt_hh
#define COOLFluiD_Physics_GReKO_NavierStokes2DGReKOPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/NavierStokesKOmegaVarSet.hh"
#include "KOmega/NavierStokesKOmegaPvt.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and Gamma-Re trasition model coupled with K-Omega turbulence model
   *
   * @author Khalil Bensassi
   */

class NavierStokes2DGReKOPuvt : public 
KOmega::NavierStokesKOmegaPvt<KOmega::NavierStokesKOmegaVarSet<NavierStokes::NavierStokes2DVarSet, 0> > {
public: // classes
  
  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DGReKOPuvt(const std::string& name,
			  Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes2DGReKOPuvt();
  
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

private :
  
  CFreal _unperturbedFluxGa;
  CFreal _unperturbedFluxRe;

}; // end of class NavierStokes2DGReKOPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokes2DGReKOPuvt_hh
