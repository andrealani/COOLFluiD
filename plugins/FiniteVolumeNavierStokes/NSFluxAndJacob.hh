#ifndef COOLFluiD_Numerics_FiniteVolume_NSFluxAndJacob_hh
#define COOLFluiD_Numerics_FiniteVolume_NSFluxAndJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NSFlux.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model
 *
 * @author Andrea Lani
 *
 */
class NSFluxAndJacob : public NSFlux<Physics::NavierStokes::NavierStokesVarSet> {
public:

  /**
   * Constructor
   */
  NSFluxAndJacob(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NSFluxAndJacob();

  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result);
  
}; // end of class NSFluxAndJacob

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NSFluxAndJacob_hh
