#ifndef COOLFluiD_Numerics_FluctSplit_STKSNCSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKSNCSchemeScalar_hh

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
class STKS_NCSchemeScalar : public STKS_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STKS_NCSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STKS_NCSchemeScalar();

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

  /// Contribution from the past to the dissipation
  std::vector<RealVector> past_diss;

/// Contribution from the past to the time integration
  std::vector<RealVector> m_time_comp;
}; // end of class STKS_NCSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKS_NCSchemeScalar_hh
