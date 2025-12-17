#ifndef COOLFluiD_Numerics_FiniteVolume_HLLDecEFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_HLLDecEFlux_hh

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
 * @author Haopeng Wang
 *
 */
class HLLDecEFlux : public FVMCC_FluxSplitter {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  HLLDecEFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HLLDecEFlux();
  
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

  /// coefficient for extra diffusion
  CFreal _kappaLax_E;

  /// use alpha coefficient
  bool _useAlpha;

  /// Add extra diffusion to HLL solver
  bool _addLax;

  /// 2D require to modify definition of plasma beta for extra diffusion
  bool _2Dornot;
  
}; // end of class HLLDecEFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HLLDecEFlux_hh
