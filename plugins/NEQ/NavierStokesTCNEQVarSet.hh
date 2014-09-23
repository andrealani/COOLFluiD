#ifndef COOLFluiD_Physics_NEQ_NavierStokesTCNEQVarSet_hh
#define COOLFluiD_Physics_NEQ_NavierStokesTCNEQVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model  for chemical NEQ
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
template <typename BASE>     
class NavierStokesTCNEQVarSet : public NavierStokesNEQVarSet<BASE> {
public: // classes
  
  typedef NavierStokesTCNEQVarSet<BASE> NEQ;
  
  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesTCNEQVarSet(const std::string& name,
			  Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesTCNEQVarSet();
  
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
  
  /**
   * Get the axisymmetric source term
   */
  virtual void getAxiSourceTerm(const RealVector& physicalData,
				const RealVector& state,
				const std::vector<RealVector*>& gradients,
				const CFreal& radius,
				RealVector& source);

  /// Compute the transport properties
  virtual void computeTransportProperties(const RealVector& state,
					  const std::vector<RealVector*>& gradients,
					  const RealVector& normal);

private:
  
  /// array to store the normal concentration gradients of species
  RealVector _normConcGradientsAxi;
  
  /// matrix to store the diffusion velocities of species multiplied by the
  /// species densities (backup)
  RealVector _rhoUdiffAxi;
  
  /// vibrational lambda
  RealVector _lambdaVibAxi;
  
}; // end of class NavierStokesTCNEQVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesTCNEQVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_NavierStokesTCNEQVarSet_hh
