#ifndef COOLFluiD_Numerics_FiniteVolume_LaxFriedFluxMultiSpecies_hh
#define COOLFluiD_Numerics_FiniteVolume_LaxFriedFluxMultiSpecies_hh

//////////////////////////////////////////////////////////////////////////////

#include "LaxFriedNSvtFlux.hh"

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
template <class UPDATEVAR>
class LaxFriedFluxMultiSpecies : public LaxFriedNSvtFlux {
public:

  /**
   * Constructor
   */
  LaxFriedFluxMultiSpecies(const std::string& name);

  /**
   * Default destructor
   */
  ~LaxFriedFluxMultiSpecies();

  /**
   * Set up private data
   */
  virtual void setup();
    
  /**
    * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);

private:

  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> _upVar;

}; // end of class LaxFriedFluxMultiSpecies

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "LaxFriedFluxMultiSpecies.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LaxFriedFluxMultiSpecies_hh
