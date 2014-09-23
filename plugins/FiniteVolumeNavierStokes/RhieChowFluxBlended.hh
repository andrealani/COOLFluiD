#ifndef COOLFluiD_Numerics_FiniteVolume_RhieChowFluxBlended_hh
#define COOLFluiD_Numerics_FiniteVolume_RhieChowFluxBlended_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/RhieChowFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the pressure-stabilized flux of Rhie-Chow corresponding
 * to the Euler physical model 2D (in conservative variables)
 *
 * @author Radek Honzatko
 *
 */
template <class UPDATEVAR>
class RhieChowFluxBlended : public RhieChowFlux<UPDATEVAR> {
public: // classes
  
  /**
   * Constructor
   */
  RhieChowFluxBlended(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RhieChowFluxBlended();
  
protected:
  
  /// Get beta speed
  virtual CFreal getBeta(const RealVector& leftData, const RealVector& rightData);
    
}; // end of class RhieChowFluxBlended

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "RhieChowFluxBlended.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RhieChowFluxBlended_hh
