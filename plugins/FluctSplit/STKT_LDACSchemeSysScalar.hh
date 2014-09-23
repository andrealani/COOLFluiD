#ifndef COOLFluiD_Numerics_FluctSplit_STKT_LDACSchemeSysScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_LDACSchemeSysScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKT_SplitterSysScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTime LDA scheme for RDS space discretization
 *
 * @author Andrea Lani
 * @author Nadege Villedieu
 *
 */
class STKT_LDACSchemeSysScalar : public STKT_SplitterSysScalar {
public:

  /**
   * Default constructor.
   */
  STKT_LDACSchemeSysScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_LDACSchemeSysScalar();

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
  
  // temporary data for computation ofthe sum of the positive part of the "space-time"upwind coeff
  
  RealMatrix m_sumKplus;
  
  RealMatrix m_invK;
  
  RealMatrix m_k;
  
  RealMatrix m_beta;
  
  RealVector m_invSumKplusScalar;
  
  RealVector m_uTemp;
  
  RealVector m_betaScalar;
  
}; // end of class SpaceTimeCsikLDASchemeSysScalar

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FluctSplit
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_LDACSchemeSysScalar_hh
