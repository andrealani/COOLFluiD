#ifndef COOLFluiD_Numerics_FluctSplit_STM_LDASchemeSysHO_hh
#define COOLFluiD_Numerics_FluctSplit_STM_LDASchemeSysHO_hh

//////////////////////////////////////////////////////////////////////////////

#include "STM_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTimeN scheme for RDS space discretization
 *
 * @author Thomas Wuilbaut
 *
 */
class STM_LDASchemeSysHO : public STM_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STM_LDASchemeSysHO(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_LDASchemeSysHO();

  /**
   * Set up
   */
  void setup();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);

  /**
   * Distribute the residual (contribution from past states)
   */ 
  void distributePast(std::vector<RealVector>& residual);

  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:
  /// storage for the sum of the kplus
  RealMatrix _sumKplus;

  /// storage for the inverse of sumKmin
  RealMatrix _invK;
 /// storage for the inverse of sumKmin
  RealMatrix m_temp_mat;

RealVector m_utemp;
}; // end of class STM_LDASchemeSysHO

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeSys_hh
