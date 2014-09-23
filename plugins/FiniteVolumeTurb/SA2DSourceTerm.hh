#ifndef COOLFluiD_Numerics_FiniteVolume_SA2DSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_SA2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "SA/NavierStokesSAVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the 2D Spalart Allmaras SourceTerm
 * variables
 *
 * @author Joao Pinto
 *
 */
class SA2DSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  SA2DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~SA2DSourceTerm();

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
  Common::SafePtr<Physics::SA::NavierStokes2DSA> _diffVarSet;
  
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
  
  /// density vector
  RealVector _rho;

  /// unperturbed Positive Part
  CFreal _unperturbedPositivePart;

  CFreal _unperturbedNegativePart;

}; // end of class SA2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SA2DSourceTerm_hh
