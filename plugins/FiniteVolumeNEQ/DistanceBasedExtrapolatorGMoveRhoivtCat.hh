#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivtCat_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivtCat_hh

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
class DistanceBasedExtrapolatorGMoveRhoivtCat : 
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
  DistanceBasedExtrapolatorGMoveRhoivtCat(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveRhoivtCat();

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
  
protected:

  /// flags indicating the nodes for which Cat has not to be applied
  std::vector<bool> _doNotApplyCat;
  
  /// names of the TRSs for which Cat has not to be applied
  std::vector<std::string> _trsWithNoCat;
  
}; // end of class DistanceBasedExtrapolatorGMoveRhoivtCat

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivtCat_hh
