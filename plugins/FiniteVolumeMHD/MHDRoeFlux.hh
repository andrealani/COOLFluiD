#ifndef COOLFluiD_Numerics_FiniteVolume_MHDRoeFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_MHDRoeFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/RoeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux corresponding to the MHD
 * physical model corresponding to the solution variable set
 * SOLUTION_VAR and to the update variable set UPDATE_VAR
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
template <class SOLUTION_VAR>
class MHDRoeFlux : public RoeFlux {
public:

  /**
   * Constructor
   */
  MHDRoeFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MHDRoeFlux();

  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
protected:

  /// acquaintance of the concrete variable set
  Common::SafePtr<SOLUTION_VAR> _solutionVarSet;
  
  /// array storing the temporary delta flux
  RealVector    _deltaFlux;
  
  /// array storing the temporary solution
  RealVector    _tempResult;
  
  /// temporary column right evectors
  RealVector   _columnRightEv;

  /// temporary wave strengths
  RealVector   _waveStrengths;

  /// rotation matrix
  RealMatrix   _rm;

  /// rotation matrix
  RealMatrix   _rmInv;

}; // end of class MHDRoeFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MHDRoeFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MHDRoeFlux_hh
