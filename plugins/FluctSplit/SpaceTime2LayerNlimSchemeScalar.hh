#ifndef COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNlimSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNlimSchemeScalar_hh

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
class SpaceTime2LayerNlimSchemeScalar : public STM_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  SpaceTime2LayerNlimSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~SpaceTime2LayerNlimSchemeScalar();

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

  /**
   * Set up
   */
  void setup();

private:

  RealVector _sumKmin;

  RealVector _sumKminU;

  RealVector _uInflow;

  RealVector _uMin;

  RealVector _uTemp;

  RealVector _phy;

  RealVector _sumBeta;

  RealVector _sumBetaLim;

  std::vector<RealVector> _resTemp;

  std::vector<RealVector> _beta;

  std::vector<RealVector> _interBeta;

  CFreal _betaLim;

  CFreal _interBetaLim;

}; // end of class SpaceTime2LayerNlimSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeScalar_hh
