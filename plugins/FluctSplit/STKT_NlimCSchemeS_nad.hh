#ifndef COOLFluiD_Numerics_FluctSplit_STKT_NlimCSchemeS_nad_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_NlimCSchemeS_nad_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKT_SplitterSys.hh"

#include "NavierStokes/EulerTerm.hh"
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTime N Limited scheme for RDS space discretization
 * We consider it as the space-Time LDA + a dissipation. This approach is developped 
 * In the thesis of Mario Ricchiuto (see pp 139)
 * @author Nadege Villedieu
 *
 */
class STKT_NlimCSchemeS_nad : public STKT_SplitterSys {
public:

  /**
   * Default constructor.
   */
  STKT_NlimCSchemeS_nad(const std::string& name);

 /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default destructor
   */
  ~STKT_NlimCSchemeS_nad();

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

  RealVector m_diss;

  RealVector m_sumKplusU;

  RealVector _phy;

  RealVector _phyChar;

  std::vector<RealVector> _residualChar;

  RealMatrix _rightEigenVector;

  RealMatrix _leftEigenVector;

  RealVector _sumBeta;

  std::vector<RealVector> _beta;

  std::vector<RealVector> _betaLim;

  RealVector m_normal;

  RealVector m_phiT;

  Common::SafePtr<Physics::NavierStokes::EulerTerm> m_cterm;
}; // end of class STKT_NlimCSchemeS_nad

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_NlimCSchemeS_nad_hh
