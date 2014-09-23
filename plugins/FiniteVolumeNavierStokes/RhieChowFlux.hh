#ifndef COOLFluiD_Numerics_FiniteVolume_RhieChowFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_RhieChowFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_FluxSplitter.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the pressure-stabilized flux of Rhie-Chow corresponding
 * to the Euler physical model (in conservative variables)
 *
 * @author Radek Honzatko
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class RhieChowFlux : public FVMCC_FluxSplitter {
public: // classes
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  RhieChowFlux(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~RhieChowFlux();
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Compute the flux in the current face
   */
  virtual void compute(RealVector& result);

protected:
  
  /// Get beta speed
  virtual CFreal getBeta(const RealVector& leftData, const RealVector& rightData);
  
protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffusiveVarSet;

  /// interface x speed component 12
  CFreal m_u12;
  
  /// interface y speed component 12
  CFreal m_v12;
  
  /// interface z speed component 12
  CFreal m_w12;
  
  /// interface density
  CFreal m_rho12;
  
  /// interface normal momentum
  CFreal m_rhoUn12;
  
  /// temporary unit normal
  RealVector m_TempUnitNormal;
  
  /// temporary average state
  RealVector _avState;
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> m_DummyGradients;
  
  /// option for pressure stabilization
  bool m_PressStab;
  
  /// sets the value of pressure dissipation scaling factor Lambda to a constant chosen value
  CFreal m_PressDissipScale;
  
  /// global flow speed estimate
  CFreal m_globalFlowSpeedEstimate;
  
}; // end of class RhieChowFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "RhieChowFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RhieChowFlux_hh
