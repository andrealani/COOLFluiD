#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMovePivt_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMovePivt_hh

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
class DistanceBasedExtrapolatorGMovePivt : public DistanceBasedExtrapolatorGMoveRhoivt {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMovePivt(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMovePivt();
  
  /**
   * Set up private data needed by the computation
   */
  virtual void setup();

protected:
  
  /**
   * Transform a given state in another one before accumulating the values
   */
  virtual void transform(const RealVector& in, RealVector& out);
  
  /**
   * Transform back a given state in the original one after accumulating the values
   */
  virtual void transformBack(RealVector& nstate);
  
  /// Set the species-related variables in a given state
  void setSpeciesVariables(RealVector& nstate);
  
protected:
  
  /// temporary partial densities
  RealVector _rhos;
  
  /// temporary species gas constants
  RealVector _Rspecies;
  
  /// R*T*y_i/M_i per each species
  RealVector _RTYovM;
  
}; // end of class DistanceBasedExtrapolatorGMovePivt

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMovePivt_hh
