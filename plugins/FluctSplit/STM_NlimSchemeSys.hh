#ifndef COOLFluiD_Numerics_FluctSplit_STM_NlimSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STM_NlimSchemeSys_hh

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

class STM_NlimSchemeSys : public STM_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STM_NlimSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_NlimSchemeSys();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);
  
  /**
   * Distribute the residual (contribution from past states)
   */
  void distributePast(const std::vector<Framework::State*>& tStates);
  
  /**
   * Set Up
   */
  void setup();
  
private:

  RealMatrix _sumKmin;

  RealMatrix _sumKplus;

  RealMatrix _invK;

  RealVector _sumKminU;

  RealVector _uInflow;

  RealMatrix _tempMat;

  RealVector _uMin;

  RealVector _uTemp;

  RealVector _phy;

  RealVector _phyChar;

  std::vector<RealVector> _residualChar;

  RealMatrix _rightEigenVector;

  RealMatrix _leftEigenVector;

  RealVector _sumBeta;

  std::vector<RealVector> _beta;

  std::vector<RealVector> _betaLim;


}; // end of class STM_NlimSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NlimSchemeSys_hh
