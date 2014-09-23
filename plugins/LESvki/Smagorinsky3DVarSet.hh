#ifndef COOLFluiD_Physics_LESvki_Smagorinsky3DVarSet_hh
#define COOLFluiD_Physics_LESvki_Smagorinsky3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LESvki {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class Smagorinsky3DVarSet : public NavierStokes::NavierStokesTurb3DVarSet {
public: // classes

  /**
   * Constructor
   * @see Smagorinsky3D
   */
   Smagorinsky3DVarSet(const std::string& name,
		       Common::SafePtr<Framework::PhysicalModelImpl> model) :
     NavierStokes::NavierStokesTurb3DVarSet(name, model)
  {
  }

  /**
   * Default destructor
   */
  virtual ~Smagorinsky3DVarSet()
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
   * Get the axisymmetric source term
   */
  virtual void getAxiSourceTerm(const RealVector& physicalData,
				const RealVector& state,
				const std::vector<RealVector*>& gradients,
				const CFreal& radius,
				RealVector& source);
  
  /**
   * Get the heat flux
   */
  virtual CFreal getHeatFlux(const RealVector& state,
                             const std::vector<RealVector*>& gradients,
                             const RealVector& normal);
  
}; // end of class Smagorinsky3DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESvki

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LESvki_Smagorinsky3DVarSet_hh
