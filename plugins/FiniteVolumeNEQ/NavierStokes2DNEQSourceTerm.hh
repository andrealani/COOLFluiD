#ifndef COOLFluiD_Numerics_FiniteVolume_NavierStokes2DNEQSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_NavierStokes2DNEQSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NavierStokes2DAxiSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 *
 */
template <class EULERVAR, class NSVAR>
class NavierStokes2DNEQSourceTerm : public NavierStokes2DAxiSourceTerm<EULERVAR, NSVAR> {
  
public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  NavierStokes2DNEQSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokes2DNEQSourceTerm();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source);
  
private:
  
  // array of gradients
  std::vector<RealVector*> m_gradients;
  
}; // end of class NavierStokes2DNEQSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DNEQSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokes2DNEQSourceTerm_hh
