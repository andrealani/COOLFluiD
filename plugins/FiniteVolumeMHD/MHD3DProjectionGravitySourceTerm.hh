#ifndef COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionGravitySourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionGravitySourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents gravitational force source term for 3D conservative
 * variables for hyperbolic divergence cleaning
 *
 * @author Mehmet Sarp Yalim
 *
 */
class MHD3DProjectionGravitySourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  MHD3DProjectionGravitySourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD3DProjectionGravitySourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
    
    _globalSockets.createSocketSink<Framework::State*>("states");
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
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// MHD physical data at the cell center
  RealVector _physicalData;

  /// x-coordinate of the center of mass of the external object if different than origin
  CFreal _xMassCenter;

  /// y-coordinate of the center of mass of the external object if different than origin
  CFreal _yMassCenter;

  /// z-coordinate of the center of mass of the external object if different than origin
  CFreal _zMassCenter;

}; // end of class MHD3DProjectionGravitySourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionGravitySourceTerm_hh
