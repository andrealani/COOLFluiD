#ifndef COOLFluiD_Numerics_FluctSplit_STKS_LDACSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STKS_LDACSchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKS_SplitterSys.hh"

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
class STKS_LDACSchemeSys : public STKS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STKS_LDACSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STKS_LDACSchemeSys();

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

  RealMatrix m_inv_K;

  RealVector m_uTemp;

  RealVector m_phiT;
  

}; // end of class STKS_LDASchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKS_LDACSchemeSys_hh
