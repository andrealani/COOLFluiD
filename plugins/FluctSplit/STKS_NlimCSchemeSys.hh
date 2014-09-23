#ifndef COOLFluiD_Numerics_FluctSplit_STKS_NlimCSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_STKS_NlimCSchemeSys_hh

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
class STKS_NlimCSchemeSys : public STKS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STKS_NlimCSchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~STKS_NlimCSchemeSys();

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

 /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);


private:


/// temporary data for computation ofthe sum of the positive part of the "space-time"upwind coeff
  RealMatrix m_sumKplus;

  RealMatrix _invK;
  
  RealVector m_uTemp;

  RealVector m_diss;

  RealVector m_sumKplusU;

  std::vector<RealVector> past_diss;

    std::vector<RealVector> m_time_comp;

  CFreal m_angle;

  RealVector m_normal;

  RealVector _phy;

  RealVector _phyChar;

  std::vector<RealVector> _residualChar;

  RealMatrix _rightEigenVector;

  RealMatrix _leftEigenVector;

  RealVector _sumBeta;

  std::vector<RealVector> _beta;

  std::vector<RealVector> _betaLim;

}; // end of class STKS_NlimCSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKS_NlimCSchemeSys_hh
