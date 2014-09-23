#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveCat_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveCat_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNEQ/DistanceBasedExtrapolatorGMoveRhoivt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object to be used 
 * in combination with BCs that moves the ghost states (nodes) like
 * @see NoSlipWallIsothermalNSPvt
 *
 * @author Andrea Lani
 *
 */
class DistanceBasedExtrapolatorGMoveCat : public DistanceBasedExtrapolatorGMoveRhoivt {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveCat(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveCat();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set up private data needed by the computation
   */
  virtual void setup();
  
protected:
    
  /**
   * Apply the boundary condition
   */
  virtual void applyBC();
    
}; // end of class DistanceBasedExtrapolatorGMoveCat

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveCat_hh
