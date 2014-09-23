#ifndef COOLFluiD_Numerics_FluctSplit_STKT_NCSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_NCSchemeScalar_hh

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
class STKT_NCSchemeScalar : public STKT_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STKT_NCSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_NCSchemeScalar();

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

  RealVector m_diss;

  RealVector m_sumKplusU;

  RealVector m_phiT;

}; // end of class STKT_NCSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_NCSchemeScalar_hh
