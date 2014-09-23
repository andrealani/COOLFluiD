#ifndef COOLFluiD_Numerics_FiniteVolume_NavierStokesKOmega3DSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_NavierStokesKOmega3DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "FiniteVolumeKOmega.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the 3D K-Omega Source Term
 *
 * @author Thomas Wuilbaut
 * @author Milan Zaloudek
 *
 */
class NavierStokesKOmega3DSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  NavierStokesKOmega3DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesKOmega3DSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);

    _sockets.createSocketSink<RealVector>("nstates");
    _sockets.createSocketSink<CFreal>("wallDistance");
  }

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
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _varSet;

  /// corresponding diffusive variable set
  Common::SafePtr<Physics::KOmega::NavierStokes3DKOmega> _diffVarSet;
  
  /// vector to store temporary result
  RealVector _temp;

  /// average State
  RealVector _avState;

  /// Euler physical data
  RealVector _physicalData;

  /// handle to the reconstructed nodal states
  Framework::DataHandle< RealVector> _nstates;

  /// handle to the wall distance
  Framework::DataHandle< CFreal> _wallDistance;

  /// array of temporary values
  RealMatrix _values;

  /// array of temporary nodal states
  std::vector<RealVector*> _states;

  /// unperturbed Positive Part
  RealVector _unperturbedPositivePart;

  /// unperturbed Negative Part
  RealVector _unperturbedNegativePart;

  ///Vector for the gradients
  std::vector<RealVector*> _gradients;

  ///temp vector for computation of gradients
  RealVector _averageFaceValue;

}; // end of class NavierStokesKOmega3DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokesKOmega3DSourceTerm_hh
