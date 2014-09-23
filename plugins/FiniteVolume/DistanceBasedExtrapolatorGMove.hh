#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMove_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMove_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DistanceBasedExtrapolator.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

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
class DistanceBasedExtrapolatorGMove : public Framework::DistanceBasedExtrapolator<CellCenterFVMData> {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMove(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMove();

  /**
   * Set up private data needed by the computation
   */
  virtual void setup();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

protected:
  
  /**
   * Apply the nodal extrapolation
   */
  virtual void extrapolate();
  
  /**
   * Check the TRS list
   */
  virtual void checkTRSList(const CFuint count);
  
  /**
   * Apply the boundary condition
   */
  virtual void applyBC();
  
protected:
  
  // arrays of flag telling of the node is on the wall
  std::valarray<bool> _isNodeToPrescribe;
  
  /// flags indicating the nodes for which radiative equilibrium has not to be applied
  std::vector<bool> _doNotApplyRadEq;
  
  /// temperature ID
  CFuint _radTID;
  
  /// flag to activate the radiative equilibrium
  bool _radEquilibrium;

  // indices of prescribed values
  std::vector<CFuint> _wValuesIdx;

  // prescribed values
  std::vector<CFreal> _wValues;
  
  /// names of the TRSs for which radiative equilibrium has not to be applied
  std::vector<std::string> _trsWithNoRadEq;
  
}; // end of class DistanceBasedExtrapolatorGMove

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMove_hh
