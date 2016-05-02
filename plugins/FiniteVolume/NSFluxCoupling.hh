#ifndef COOLFluiD_Numerics_FiniteVolume_NSFluxCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_NSFluxCoupling_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/NSFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model with weak coupling
 *
 * @author Andrea Lani
 *
 */
template <typename DIFFVS>      
class NSFluxCoupling : public NSFlux<DIFFVS> {
public:

  /**
   * Constructor
   */
  NSFluxCoupling(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NSFluxCoupling();
  
  /**
    * Compute the flux in the current face
    */
  virtual void computeFlux(RealVector& result);
  
}; // end of class NSFluxCoupling
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/NSFluxCoupling.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NSFluxCoupling_hh
