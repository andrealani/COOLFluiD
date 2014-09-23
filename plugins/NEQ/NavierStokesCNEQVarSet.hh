#ifndef COOLFluiD_Physics_NEQ_NavierStokesCNEQVarSet_hh
#define COOLFluiD_Physics_NEQ_NavierStokesCNEQVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model for chemical NEQ
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
template <typename BASE>     
class NavierStokesCNEQVarSet : public NavierStokesNEQVarSet<BASE> {
public: // classes
  
  typedef NavierStokesCNEQVarSet<BASE> NEQ;
  
  /**
   * Constructor
   */
  NavierStokesCNEQVarSet(const std::string& name,
			 Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesCNEQVarSet();
  
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

  /// Compute the transport properties
  virtual void computeTransportProperties(const RealVector& state,
					  const std::vector<RealVector*>& gradients,
					  const RealVector& normal);
  
}; // end of class NavierStokesCNEQVarSet
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesCNEQVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_NavierStokesCNEQVarSet_hh
