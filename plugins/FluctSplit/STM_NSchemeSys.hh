#ifndef COOLFluiD_Numerics_FluctSplit_STM_NSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STM_NSchemeSys_hh

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
class STM_NSchemeSys : public STM_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STM_NSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_NSchemeSys();

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
  
private:

  RealMatrix _sumKmin;

  RealMatrix _sumKplus;

  RealMatrix _invK;

  RealVector _sumKminU;

  RealVector _uInflow;

  RealMatrix _tempMat;

  RealVector _uMin;

  RealVector _uTemp;

}; // end of class STM_NSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeSys_hh
