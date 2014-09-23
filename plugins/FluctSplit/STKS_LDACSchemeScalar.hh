#ifndef COOLFluiD_Numerics_FluctSplit_STKS_LDACSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STKS_LDACSchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKS_SplitterScalar.hh"

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
class STKS_LDACSchemeScalar : public STKS_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STKS_LDACSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STKS_LDACSchemeScalar();

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

  RealVector m_phiT;
  
}; // end of class SpaceTimeCaLDASchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_LDACSchemeScalar_hh
