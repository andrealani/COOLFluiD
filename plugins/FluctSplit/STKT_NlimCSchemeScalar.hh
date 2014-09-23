#ifndef COOLFluiD_Numerics_FluctSplit_STKT_NlimCSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_NlimCSchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKT_SplitterScalar.hh"

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
class STKT_NlimCSchemeScalar : public STKT_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STKT_NlimCSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_NlimCSchemeScalar();

  /**
   * Set up
   */
  void setup();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);

 
  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

/// temporary data for computation ofthe sum of the positive part of the "space-time"upwind coeff
  RealVector m_sumKplus;

  RealVector m_uTemp;

  RealVector m_sumKplusU;

  RealVector m_diss;

  RealVector m_phiT;
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
}; // end of class STKT_NlimCSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_NlimCSchemeScalar_hh
