#ifndef COOLFluiD_Numerics_FluctSplit_STKT_BCSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_BCSchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/STKT_SplitterSys.hh"

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
class STKT_BCSchemeSys : public STKT_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STKT_BCSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_BCSchemeSys();

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
  RealMatrix m_sumKplus;

  RealMatrix _invK;

  RealVector m_uTemp;

  std::vector<RealVector> m_diss;

  RealVector m_sumKplusU;

  std::vector<RealVector> m_phiN;

  RealVector m_phitot;

  RealVector m_absphitot;

  RealVector m_phiT;

}; // end of class STKT_BCSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_BCSchemeSys_hh
