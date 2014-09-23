#ifndef COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionJupiterSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionJupiterSourceTerm_hh

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
 * This class represents a source term for including the effect of the ionosphere
 * into the Jovian magnetosphere for projection scheme
 * for 3D conservative variables
 *
 * @author Mehmet Sarp Yalim
 *
 */
class MHD3DProjectionJupiterSourceTerm : public ComputeSourceTermFVMCC {

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
  MHD3DProjectionJupiterSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~MHD3DProjectionJupiterSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);

    _sockets.createSocketSink<Framework::State*, Framework::GLOBAL>("states");
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Set the variable set
   * @pre the input pointer is non const to allow dynamic_cast
   */
  void setVarSet(Common::SafePtr<Framework::ConvectiveVarSet> varSet);

  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
         RealVector& source);

private: // data

  /// corresponding variable set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  RealVector _physicalData;
  
  /// vector of state coordinates
  RealVector _stateCoord;

  /// distance of center of the torus from the origin
  CFreal _a0;

  /// radius of the torus
  CFreal _radiusTorus;

  /// corotation velocity of the torus (tangential)
  CFreal _velTangentialTorus;

  /// mass density production rate in the torus
  CFreal _dRhodtTorusPlasma;

  /// Temperature of the neutrals in the torus
  CFreal _torusTemp;

  /// Atomic mass of the ions in the torus
  CFreal _torusAtomicMass;
  
  /// Rotation period of the planet
  CFreal _ionoRotPeriod;

  /// Maximum Ion-Neutral collision frequency in the ionosphere
  CFreal _ionoCollFreq;

  /// Decay constant (it decreases exponentially) for the ionosphere
  CFreal _ionoDecay;

  /// Atomic mass of the ions in the ionosphere
  CFreal _ionoAtomicMass;
  
  /// Radius of the inner boundary
  CFreal _radiusBnd;

  /// Mi/Mn - Mass ratio in the ionosphere
  CFreal _ionoIonNeutralMassRatio;

  /// Temperature of the neutrals in the ionosphere
  CFreal _ionoTemp;

  /// G.M - Gravitational constant times the mass of the planet
  CFreal _gravityTimesJupiMass;
  

}; // end of class MHD3DProjectionJupiterSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHD3DProjectionJupiterSourceTerm_hh
