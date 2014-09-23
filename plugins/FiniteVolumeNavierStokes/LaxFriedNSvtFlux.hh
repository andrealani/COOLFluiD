#ifndef COOLFluiD_Numerics_FiniteVolume_LaxFriedNSvtFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_LaxFriedNSvtFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/LaxFriedFlux.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Lax-Friedrichs flux
 *
 * @author Andrea Lani
 *
 */
class LaxFriedNSvtFlux : public LaxFriedFlux {
public:

  /**
   * Constructor
   */
  LaxFriedNSvtFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LaxFriedNSvtFlux();

  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

protected:
  
  /**
   * Compute the artificial diffusion reduction coefficient
   */
  CFreal getReductionCoeff();
  
private:
  
  /// pointer to the Navier Stokes difffusive VarSet
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_nsVarSet;
  
  /// dummy gradients
  std::vector<RealVector*> m_dummyGradients;
  
  /// average state
  RealVector m_avState;
  
  /// minimum value for cell Reynolds number
  CFreal m_reynoldsMin;
  
  /// array storing the IDs of the variables corresponding to velocity
  /// components
  std::vector<CFuint> m_velocityIDs;
  
}; // end of class LaxFriedNSvtFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LaxFriedNSvtFlux_hh
