#ifndef COOLFluiD_Numerics_FiniteVolume_LaxFriedCouplingFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_LaxFriedCouplingFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Lax-Friedrichs flux
 *
 * @author Andrea Lani
 * @author Sarp Yalim
 *
 */
class LaxFriedCouplingFlux : public FVMCC_FluxSplitter {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  LaxFriedCouplingFlux(const std::string& name);

  /**
   * Default destructor
   */
  ~LaxFriedCouplingFlux();

  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
private:

  /**
   * Compute the flux in the perturbed case
   */
  void computePerturbCase(RealVector& result);

  /**
   * Compute the flux in the unperturbed case
   */
  void computeUnperturbCase(RealVector& result);

private:

  /// array storing the sum of the right and left flux
  RealVector   _sumFlux;

  /// array storing the temporary eigenvalues
  RealVector    _aDiffVec;

  /// array storing the temporary right eigenvalues
  RealVector    _rightEv;

  /// array storing the temporary left eigenvalues
  RealVector    _leftEv;

  /// temporary unit normal
  RealVector    _tempUnitNormal;

  /// diffusion reduction coefficient defined interactively
  CFreal            _currentDiffRedCoeff;

}; // end of class LaxFriedCouplingFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LaxFriedCouplingFlux_hh
