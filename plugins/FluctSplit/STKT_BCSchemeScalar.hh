#ifndef COOLFluiD_Numerics_FluctSplit_STKT_BCSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_BCSchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/STKT_SplitterScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTime B scheme for RDS space discretization
 * We consider it as the Space-Time LDA + theta * dissipation used for N scheme.
 * This approach is developped
 * In the thesis of Mario Ricchiuto (see pp 139)
 * @author Nadege Villedieu
 *
 */
class STKT_BCSchemeScalar : public STKT_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STKT_BCSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_BCSchemeScalar();

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

  std::vector<RealVector> m_diss;

  std::vector<RealVector> m_phiN;

  RealVector m_phitot;

  RealVector m_absphitot;

  RealVector m_phiT;

}; // end of class STKT_BCSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_BCSchemeScalar_hh
