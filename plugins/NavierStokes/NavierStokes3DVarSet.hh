#ifndef COOLFluiD_Physics_NavierStokes_NavierStokes3DVarSet_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokes3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class NavierStokes3DVarSet : public NavierStokesVarSet {
public: // classes

  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
    NavierStokesVarSet(name, model)
  {
  }

  /**
   * Default destructor
   */
  virtual ~NavierStokes3DVarSet()
  {
  }

  /**
   * Set up private data
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
  virtual RealMatrix& getFlux(const RealVector& state,
			      const std::vector<RealVector*>& gradients,
			      const CFreal& radius);
  
}; // end of class NavierStokes3DVarSet
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokes3DVarSet_hh
