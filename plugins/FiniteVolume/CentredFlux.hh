#ifndef COOLFluiD_Numerics_FiniteVolume_CentredFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_CentredFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes a pure central flux
 *
 * @author Kris Van den Abeele
 * @author Ghader Ghorbaniasl
 *
 */
class CentredFlux : public FVMCC_FluxSplitter {
public:
	
  /**
   * Constructor
   */
  CentredFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CentredFlux();

  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
protected:
  
  /// array storing the sum of the right and left flux
  RealVector   _sumFlux;
  
  /// temporary unit normal
  RealVector    _tempUnitNormal;
  
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
}; // end of class CentredFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CentredFlux_hh
