#ifndef COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNlimSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_SpaceTime2LayerNlimSchemeSys_hh

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
class SpaceTime2LayerNlimSchemeSys : public STM_SplitterSys {
public:

  /**
   * Default constructor.
   */
  SpaceTime2LayerNlimSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~SpaceTime2LayerNlimSchemeSys();

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
   * Sets up the scheme
   */
  void setup();

private:

  RealMatrix _sumKmin;

  RealMatrix _invK;

  RealVector _sumKminU;

  RealVector _uInflow;

  RealMatrix _tempMat;

  RealVector _uMin;

  RealVector _uTemp;

  RealVector _phy;

  RealVector _phyChar;

  std::vector<RealVector> _residualChar;

  std::vector<RealVector> _interResidualChar;

  RealMatrix _rightEigenVector;

  RealMatrix _leftEigenVector;

  RealVector _sumBeta;

  RealVector _sumBetaLim;

  std::vector<RealVector> _resTemp;

  std::vector<RealVector> _resTemp2;

  std::vector<RealVector> _beta;

  std::vector<RealVector> _interBeta;

  CFreal _betaLim;

  CFreal _interBetaLim;

}; // end of class SpaceTime2LayerNlimSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeSys_hh
