#ifndef COOLFluiD_Numerics_FiniteVolume_HLLETanakaFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_HLLETanakaFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/HLLEFlux.hh"

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
class HLLETanakaFlux : public HLLEFlux<UPDATEVAR> {
public:

  /**
   * Constructor
   */
  HLLETanakaFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HLLETanakaFlux();
  
  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
    
}; // end of class HLLETanakaFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "HLLETanakaFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HLLETanakaFlux_hh
