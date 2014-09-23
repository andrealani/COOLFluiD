#ifndef COOLFluiD_Numerics_FiniteVolume_HLLEFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_HLLEFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/HLLFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the HLLE flux corresponding to the
 * physical model defined in the template parameter VARSET
 *
 * @author Andrea Lani
 *
 */
template <typename UPDATEVAR>      
class HLLEFlux : public HLLFlux {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  HLLEFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HLLEFlux();
  
  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
protected:
   
  /// update variable set
  Common::SafePtr<UPDATEVAR> _upVar;
  
  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStates;
  
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
  /// use Roe average
  CFuint _useRoeAverage;
  
}; // end of class HLLEFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "HLLEFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HLLEFlux_hh
