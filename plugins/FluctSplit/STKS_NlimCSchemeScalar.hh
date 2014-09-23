#ifndef COOLFluiD_Numerics_FluctSplit_STKS_NlimCSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKS_NlimCSchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKS_SplitterScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTime N scheme for RDS space discretization
 * We consider it as the Pace-Time LDA + a dissipation. This approach is developped 
 * In the thesis of Mario Ricchiuto (see pp 139)
 * @author Nadege Villedieu
 *
 */
class STKS_NlimCSchemeScalar : public STKS_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STKS_NlimCSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STKS_NlimCSchemeScalar();

  /**
   * Set up
   */
  void setup();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);

 /**
*
* Compute dissipation from the past
*/
  void ComputePastDissipationAndTimeComp( const std::vector<Framework::State*>& tStates );


  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:


/// temporary data for computation ofthe sum of the positive part of the "space-time"upwind coeff
  RealVector m_sumKplus;

  RealVector m_uTemp;

  RealVector m_diss;

  RealVector m_sumKplusU;

  std::vector<RealVector> past_diss;

    std::vector<RealVector> m_time_comp;

  /// the residual of the standard N scheme before limiting is applied
  std::vector<RealVector> m_phiN;

  /// the beta coefficient relative to the unlimited N scheme
  std::vector<RealVector> m_betaN;

  /// the positive values of the beta coefficients if the unlimited N scheme
  std::vector<RealVector> m_betaNPlus;

  /// residual total of the element
  RealVector m_phitot;
  
  /// the sum of the m_betaPlus over all the states
  RealVector m_sumBetaNPlus;

  /// the final distribution coefficient for the limited N scheme
  std::vector<RealVector> m_beta;
}; // end of class STKS_NlimCSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKS_NlimCSchemeScalar_hh
