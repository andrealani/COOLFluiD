#ifndef COOLFluiD_Numerics_FiniteVolume_PoissonNEQSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_PoissonNEQSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BaseDataSocketSink;
  }
  
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
class PoissonNEQSourceTerm : public ComputeSourceTermFVMCC {
  
public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  PoissonNEQSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~PoissonNEQSourceTerm();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private:

  /// corresponding variable set
  Common::SafePtr<EULERVAR> m_varSet;
  
  /// corresponding diffusive variable set
  Common::SafePtr<NSVAR> m_diffVarSet;
  
  /// socket for B field (Bx, By, Bz)
  Framework::DataSocketSink<CFreal> socket_Bfield;
  
  /// Euler physical data
  RealVector m_physicalData;
  
}; // end of class PoissonNEQSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "PoissonNEQSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PoissonNEQSourceTerm_hh
