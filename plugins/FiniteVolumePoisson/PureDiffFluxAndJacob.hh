#ifndef COOLFluiD_Numerics_FiniteVolume_PureDiffFluxAndJacob_hh
#define COOLFluiD_Numerics_FiniteVolume_PureDiffFluxAndJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumePoisson/PureDiffFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux and analytical flux jacobian 
 * corresponding to the Poisson model
 *
 * @author Andrea Lani
 *
 */
class PureDiffFluxAndJacob : public PureDiffFlux {
public:
  
  /**
   * Constructor
   */
  PureDiffFluxAndJacob(const std::string& name);

  /**
   * Default destructor
   */
  ~PureDiffFluxAndJacob();
  
  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result);

    
}; // end of class PureDiffFluxAndJacob
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PureDiffFluxAndJacob_hh
