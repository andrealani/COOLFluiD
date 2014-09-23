#ifndef COOLFluiD_Numerics_FiniteVolume_LaxFriedBCCorrFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_LaxFriedBCCorrFlux_hh

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
 *
 */
class LaxFriedBCCorrFlux : public FVMCC_FluxSplitter {
public:

  /**
   * Constructor
   */
  LaxFriedBCCorrFlux(const std::string& name);

  /**
   * Default destructor
   */
  ~LaxFriedBCCorrFlux();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Unset up private data
   */
  void unsetup();

  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
private:

  /// array storing the sum of the right and left flux
  RealVector   _sumFlux;

  /// array storing the temporary right eigenvalues
  RealVector    _rightEv;

  /// array storing the temporary left eigenvalues
  RealVector    _leftEv;

  /// temporary unit normal
  RealVector    _tempUnitNormal;

  /// average State
  Framework::State* _avState;

  /// physical data  array for the average State
  RealVector _avPdata;

  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;

}; // end of class LaxFriedBCCorrFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LaxFriedBCCorrFlux_hh
