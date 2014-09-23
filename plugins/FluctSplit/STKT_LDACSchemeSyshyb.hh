#ifndef COOLFluiD_Numerics_FluctSplit_STKT_LDACSchemeSyshyb_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_LDACSchemeSyshyb_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKT_SplitterSyshyb.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTime LDA scheme for RDS space discretization
 *
 * @author Nadege Villedieu
 *
 */
class STKT_LDACSchemeSyshyb : public STKT_SplitterSyshyb {
public:

  /**
   * Default constructor.
   */
  STKT_LDACSchemeSyshyb(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_LDACSchemeSyshyb();

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

  RealVector m_phiT;

  RealMatrix _k;

  RealMatrix _beta;

}; // end of class SpaceTimeCsikLDASchemeSyshyb

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_LDACSchemeSyshyb_hh
