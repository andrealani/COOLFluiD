#ifndef COOLFluiD_Numerics_FiniteVolume_Rotation3DSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_Rotation3DSourceTerm_hh

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
      class Euler3DRotationVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Rotation source term for 3D Navier-Stokes
 *
 * @author Tom Vestraete
 *
 */
class Rotation3DSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  Rotation3DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~Rotation3DSourceTerm();

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
  Common::SafePtr<Physics::NavierStokes::Euler3DRotationVarSet> _varSet;

  /// vector to store temporary result
  RealVector _temp;

  /// Euler physical data
  RealVector _physicalData;


}; // end of class Rotation3DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Rotation3DSourceTerm_hh
