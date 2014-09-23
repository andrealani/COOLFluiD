#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivtLTE_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivtLTE_hh

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
 * @see NoSlipWallIsothermalNSvt
 *
 * @author Andrea Lani
 *
 */
class DistanceBasedExtrapolatorGMoveRhoivtLTE : 
	public DistanceBasedExtrapolatorGMoveRhoivt {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveRhoivtLTE(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveRhoivtLTE();

  /**
   * Set up private data needed by the computation
   */
  virtual void setup();

protected:
  
  /**
   * Apply the nodal extrapolation
   */
  virtual void extrapolate();
  
  /**
   * Transform back a given state in the original one after accumulating the values
   */
  virtual void transformBack(RealVector& nstate);
  
  /**
   * Apply the boundary condition
   */
  virtual void applyBC();
  
protected:

  /// array of species molar fractions
  RealVector _xs;
  
  /// flags indicating the nodes for which LTE has not to be applied
  std::vector<bool> _doNotApplyLTE;
  
  /// names of the TRSs for which LTE has not to be applied
  std::vector<std::string> _trsWithNoLTE;
  
}; // end of class DistanceBasedExtrapolatorGMoveRhoivtLTE

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivtLTE_hh
