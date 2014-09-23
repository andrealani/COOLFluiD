#ifndef COOLFluiD_Numerics_FiniteVolume_HLLFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_HLLFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the HLL flux corresponding to the
 * physical model defined in the template parameter VARSET
 *
 * @author Michel Rasquin
 *
 */
class HLLFlux : public FVMCC_FluxSplitter {
public:

  /**
   * Constructor
   */
  HLLFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HLLFlux();
  
  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
protected:
  
  /// array storing the sum of the right flux
  RealVector   _rightFlux;

 /// array storing the sum of the left flux
  RealVector   _leftFlux;

  /// array storing the temporary right eigenvalues
  RealVector    _rightEv;

  /// array storing the temporary left eigenvalues
  RealVector    _leftEv;

  /// temporarry unit normal
  RealVector    _tempUnitNormal;
  
  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStates;
  
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;

}; // end of class HLLFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HLLFlux_hh
