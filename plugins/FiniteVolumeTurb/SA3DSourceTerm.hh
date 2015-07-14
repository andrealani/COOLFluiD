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
 * This class represents the 3D Spalart Allmaras - noft2 SourceTerm and as a user option the compressibility correction term SourceTerm
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
  virtual ~SA3DSourceTerm();

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
  virtual void setup();

  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
  
protected:
  
  /// compute the density spatial derivatives when LS gradients cannot be used
  virtual void computeDRhoDX(Framework::GeometricEntity *const element);
  
protected: //methods
  
  ///@return the distance to the wall
  virtual CFreal getDistance (Framework::GeometricEntity *const element);

  // this method is needed for the SA - Comp model and for the DDES and IDDES modes
  ///@return the sum of the velocity gradients
  virtual CFreal compSumOfVelocityGrads ();
  
  //CG: this method is added so as for the DDES and IDDES modes to have access
  ///@return the laminar kinematic viscosity
  virtual CFreal getLamViscosity()
  {
    return _NIU;
  };
  
  //CG: this method ia added so as for the DDES and IDDES modes to have access
  ///@return the turbulent kinematic viscosity
  virtual CFreal getTurbViscosity()
  {
    return _NIUturbulent;
  };
  
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _varSet;
  
  /// corresponding diffusive variable set
  Common::SafePtr<Physics::SA::NavierStokes3DSA> _diffVarSet;
  
  /// vector to store temporary result
  RealVector _avStateFace;

  /// Euler physical data
  RealVector _physicalData;
  
  /// Euler physical data
  RealVector _physicalDataFace;
  
  /// handle to the reconstructed nodal states
  Framework::DataHandle< RealVector> _nstates;

  /// handle to the wall distance
  Framework::DataHandle< CFreal> _wallDistance;
  
  /// unperturbed Positive Part
  CFreal _unperturbedPositivePart;
  
  CFreal _unperturbedNegativePart;
  
private:
  
  // turbulent and laminar viscosity
  CFreal _NIUtilda;
  CFreal _NIU;
  CFreal _NIUturbulent;
  
  // distance to the wall
  CFreal _d;
  
  // gradients of density 
  CFreal _dRhodX;
  CFreal _dRhodY;
  CFreal _dRhodZ; 
  
  // gradients of velocities
  CFreal _dVdX ;
  CFreal _dUdY ;
  CFreal _dUdZ ; 
  CFreal _dWdX ; 
  CFreal _dVdZ ; 
  CFreal _dWdY ; 
  
  // these terms are added for the comp model
  CFreal _dUdX ; 
  CFreal _dVdY ;
  CFreal _dWdZ;
  
  // this variable is flag for the SA - comp model
  bool _CompTerm;
  
  // tell if a perface un-reactive gas is assumed
  bool _isPerfectGas;
  
}; // end of class SA3DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SA3DSourceTerm_hh
