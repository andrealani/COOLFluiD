#ifndef COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNSchemeSys_hh

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
class SpaceTime2LayerNSchemeSys : public STM_SplitterSys {
public:

  /**
   * Default constructor.
   */
  SpaceTime2LayerNSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~SpaceTime2LayerNSchemeSys();

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
  void distributeInterK1(const std::vector<Framework::State*>& tStates,
                                        std::vector<RealVector>& residual);

  /**
   * Distribute the residual
   */
  void distributeInterK2(const std::vector<Framework::State*>& tStates,
                                        std::vector<RealVector>& residual);

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);
  
private:

  RealMatrix _sumKmin;

  RealMatrix _invK;

  RealVector _sumKminU;

  RealVector _uInflow;

  RealMatrix _tempMat;

  RealVector _uMin;

  RealVector _uTemp;

}; // end of class SpaceTime2LayerNSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeSys_hh
