#ifndef COOLFluiD_Numerics_FluctSplit_STM_LDACSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STM_LDACSchemeSys_hh

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
class STM_LDACSchemeSys : public STM_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STM_LDACSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_LDACSchemeSys();

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
  void distributePast(const std::vector<Framework::State*>& tStates);

  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealMatrix _sumKmin;

  RealMatrix _sumKplus;

  RealVector _uMin;
  RealVector _temp;
  RealMatrix _uTemp;
  RealMatrix _temp_mat;

  /// storage for the inverse of sumKmin
  RealMatrix _invK;

  /// storage for the identity matrix (only the diagonal is stored)
  RealVector _identity;

}; // end of class STM_LDACSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDACSchemeSys_hh
