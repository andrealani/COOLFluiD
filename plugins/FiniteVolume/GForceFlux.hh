#ifndef COOLFluiD_Numerics_FiniteVolume_GForceFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_GForceFlux_hh

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
class GForceFlux : public FVMCC_FluxSplitter {
public:
  
  /**
   * Constructor
   */
  GForceFlux(const std::string& name);

  /**
   * Default destructor
   */
  ~GForceFlux();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Set up private data
   */
  virtual void unsetup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
private:
  
  /**
   * Compute the blending coefficient corresponding to the given state
   */
  CFreal computeBlendingCoeff(Framework::State& state);
  
private:
  
  /// array storing the right flux
  RealVector   _fluxR;

  /// array storing the left flux
  RealVector   _fluxL;

  /// array storing the temporary right eigenvalues
  RealVector    _rightEv;

  /// array storing the temporary left eigenvalues
  RealVector    _leftEv;

  /// temporary unit normal
  RealVector    _tempUnitNormal;

  /// gradient
  RealVector    _gradient;

  /// blending coefficient array
  RealVector    _bCoeffArray;

  /// state for the lax wendroff scheme
  Framework::State* _stateLW;

  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
  /// physical data array
  RealVector _pdataLW;
  
  /// blending coefficient defined interactively
  CFreal _coeff;

  /// maximum blending reduction coefficient defined interactively
  CFreal _maxCoeff;

  /// gradient variable ID for computing the diffusion coefficient
  CFuint _gradVarID;

  /// threshold residual to control the freezing of the blending
  /// coefficient
  CFreal _freezeCoeffRes;

}; // end of class GForceFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_GForceFlux_hh
