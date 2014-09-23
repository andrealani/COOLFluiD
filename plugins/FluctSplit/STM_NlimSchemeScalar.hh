#ifndef COOLFluiD_Numerics_FluctSplit_STM_NlimSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STM_NlimSchemeScalar_hh

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
class STM_NlimSchemeScalar : public STM_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STM_NlimSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_NlimSchemeScalar();

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

  RealVector _sumKmin;

  RealVector _sumKminU;

  RealVector _sumKplus;

  RealVector _uMin;

  RealVector _uTemp;

  RealVector _uInflow;

  RealVector _phy;

  RealVector _sumBeta;

  RealVector _sumBetaLim;

  std::vector<RealVector> _beta;

  std::vector<RealVector> _betaLim;

}; // end of class STM_NlimSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NlimSchemeScalar_hh
