#ifndef COOLFluiD_Physics_NavierStokes_NavierStokes2DVarSet_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokes2DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class NavierStokes2DVarSet : public NavierStokesVarSet {
public: // classes

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DVarSet(const std::string& name, 
		       Common::SafePtr<Framework::PhysicalModelImpl> model) :
    NavierStokesVarSet(name, model)
  {
  }
  
  /**
   * Default destructor
   */
  virtual ~NavierStokes2DVarSet()
  {
  }

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);
  
  /**
   * Get the axisymmetric source term
   */
  virtual void getAxiSourceTerm(const RealVector& physicalData,
				const RealVector& state,
				const std::vector<RealVector*>& gradients,
				const CFreal& radius,
				RealVector& source);
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);
  
}; // end of class NavierStokes2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokes2DVarSet_hh
