#ifndef COOLFluiD_Numerics_FluctSplit_STM_LDASchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STM_LDASchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "STM_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTimeLDA scheme for RDS space discretization
 *
 * @author Thomas Wuilbaut
 *
 */
class STM_LDASchemeSys : public STM_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STM_LDASchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_LDASchemeSys();

  /**
   * Set up
   */
  void setup();

  /**
   * Distribute the residual (contribution from past states)
   */
  void distributePast(const std::vector<Framework::State*>& tStates);

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);                  

private:

  /// storage for the sum of Kmin
  RealMatrix _sumKmin;

  RealMatrix _sumKplus;

  /// storage for the inverse of sumKmin
  RealMatrix _invK;

  /// storage for temporary vector
  RealVector _uMin;

  /// storage for temporary matrix
  RealMatrix _tempMat;

  /// storage for temporary vector
  std::vector<RealVector> _phyt;

  /// storage for temporary vector
  RealVector _phys;

  /// storage for the identity matrix (only the diagonal is stored)
  RealVector _identity;

}; // end of class STM_LDASchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeSys_hh
