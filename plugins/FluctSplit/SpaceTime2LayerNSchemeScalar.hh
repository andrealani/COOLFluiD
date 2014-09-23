#ifndef COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNSchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "STM_SplitterScalar.hh"

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
class SpaceTime2LayerNSchemeScalar : public STM_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  SpaceTime2LayerNSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~SpaceTime2LayerNSchemeScalar();

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

  RealVector _sumKmin;

  RealVector _sumKminU;

  RealVector _uInflow;

  RealVector _uMin;

}; // end of class SpaceTime2LayerNSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeScalar_hh
