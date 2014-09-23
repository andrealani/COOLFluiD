#ifndef COOLFluiD_Physics_NavierStokes_NavierStokes1DVarSet_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokes1DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 1D for primitive
   * variables
   *
   * @author Andrea Lani
   * @author Alessandro Munafo'
   */
class NavierStokes1DVarSet : public NavierStokesVarSet {
public: // classes

  /**
   * Constructor
   * @see NavierStokes1D
   */
  NavierStokes1DVarSet(const std::string& name,
		       Common::SafePtr<Framework::PhysicalModelImpl> model) :
    NavierStokesVarSet(name, model)
  {
  }

  /**
   * Default destructor
   */
  virtual ~NavierStokes1DVarSet()
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
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);
  
  /**
   * Get the heat flux
   */
  virtual CFreal getHeatFlux(const RealVector& state,
                             const std::vector<RealVector*>& gradients,
                             const RealVector& normal);
  
}; // end of class NavierStokes1DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokes1DVarSet_hh
