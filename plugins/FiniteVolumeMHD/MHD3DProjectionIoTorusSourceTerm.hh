#ifndef COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionIoTorusSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionIoTorusSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

// #include "Environment/FileHandlerOutput.hh"
#include "Common/SafePtr.hh"

// #include "Environment/SingleBehaviorFactory.hh"
#include "Framework/State.hh"

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class GeometricEntity; }

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term for including the effect of the torus of Io
 * into the Jovian magnetosphere for projection scheme
 * for 3D conservative variables
 *
 * @author Mehmet Sarp Yalim
 *
 */
class MHD3DProjectionIoTorusSourceTerm : public ComputeSourceTermFVMCC {

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  MHD3DProjectionIoTorusSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD3DProjectionIoTorusSourceTerm();

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

  /// vector of state coordinates
  RealVector _stateCoord;

  /// distance of center of the torus from the origin
  CFreal _a0;

  /// radius of the torus
  CFreal _radiusTorus;

  /// corotation velocity of the torus (tangential)
  CFreal _velTangentialTorus;

  /// the rate of change of plasma density in the torus
  CFreal _dRhodtTorusPlasma;

}; // end of class MHD3DProjectionIoTorusSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionIoTorusSourceTerm_hh
