#ifndef COOLFluiD_Numerics_FiniteVolume_SA3DSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_SA3DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "SA/NavierStokesSAVarSetTypes.hh" // incudes the NSVarSet and NSPuvt

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
 * This class represents the 3D Spalart Allmaras - noft2 SourceTerm
 * for more details see Allmaras, S. R., Johnson, F. T., and Spalart, P. R., "Modifications and Clarifications for the Implementation 
 * of the Spalart-Allmaras Turbulence Model," ICCFD7-1902, 7th International Conference on Computational Fluid Dynamics, Big Island, Hawaii, 9-13 July 2012.	
 *
 * @author Joao Pinto
 * @modified Christos Gkoudesnes
 *
 */
class SA3DSourceTerm : public ComputeSourceTermFVMCC {

public:
  
  /**
   * Defines the Config Option's of this class
   * @param options an OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  SA3DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~SA3DSourceTerm();

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
  Common::SafePtr<Physics::SA::NavierStokes3DSA> _diffVarSet;
  
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
  
  /// unperturbed Positive Part
  CFreal _unperturbedPositivePart;

  CFreal _unperturbedNegativePart;
  
  // this variable is flag for the SA - comp model
  bool _CompTerm;

}; // end of class SA3DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SA3DSourceTerm_hh
