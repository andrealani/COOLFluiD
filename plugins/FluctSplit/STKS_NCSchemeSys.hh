#ifndef COOLFluiD_Numerics_FluctSplit_STKS_NCSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STKS_NCSchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKS_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTime N scheme for RDS space discretization
 * We consider it as the space-Time LDA + a dissipation. This approach is developped 
 * In the thesis of Mario Ricchiuto (see pp 139)
 * @author Nadege Villedieu
 *
 */
class STKS_NCSchemeSys : public STKS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STKS_NCSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STKS_NCSchemeSys();

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
/**
*
* Compute dissipation from the past
*/
  void ComputePastDissipationAndTimeComp( const std::vector<Framework::State*>& tStates );




private:


/// temporary data for computation ofthe sum of the positive part of the "space-time"upwind coeff
  RealMatrix m_sumKplus;
  
  RealVector m_sumKplusU;

  RealVector m_uTemp;

  RealVector m_diss;

  RealMatrix _invK;
  
  std::vector<RealVector> past_diss;

    std::vector<RealVector> m_time_comp;
}; // end of class STKS_NCSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKS_NCSchemeSys_hh
