#ifndef COOLFluiD_Numerics_FluctSplit_STKT_LDACSchemeScalarQuad_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_LDACSchemeScalarQuad_hh

//////////////////////////////////////////////////////////////////////////////

#include "STKT_SplitterScalarQuad.hh"

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
class STKT_LDACSchemeScalarQuad : public STKT_SplitterScalarQuad {
public:

  /**
   * Default constructor.
   */
  STKT_LDACSchemeScalarQuad(const std::string& name);

  /**
   * Default destructor
   */
  ~STKT_LDACSchemeScalarQuad();

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
  
}; // end of class SpaceTimeCsikLDASchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics__STKT_LDACSchemeScalarQuad_hh
