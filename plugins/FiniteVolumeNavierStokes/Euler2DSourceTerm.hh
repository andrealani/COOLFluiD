#ifndef COOLFluiD_Numerics_FiniteVolume_Euler2DSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_Euler2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
    }
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
class Euler2DSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  Euler2DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DSourceTerm();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
    
private: // data
  /// corresponding variable set
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> _varSet;

  /// vector to store temporary result
  RealVector _temp;

  /// Euler physical data
  RealVector _physicalData;

}; // end of class Euler2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Euler2DSourceTerm_hh
