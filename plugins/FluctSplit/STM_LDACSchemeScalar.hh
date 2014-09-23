#ifndef COOLFluiD_Numerics_FluctSplit_STM_LDACSchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STM_LDACSchemeScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "STM_SplitterScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SpaceTimeN scheme for RDS space discretization
 *
 * @author Thomas Wuilbaut
 *
 */
class STM_LDACSchemeScalar : public STM_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STM_LDACSchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_LDACSchemeScalar();

  /**
   * Set up
   */
  void setup();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);

  /**
   * Distribute the residual (contribution from past states)
   */
  void distributePast(const std::vector<Framework::State*>& tStates);

  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealVector _sumKmin;

  RealVector _sumKplus;

  RealVector _uMin;

  RealVector _uTemp;

}; // end of class STM_LDACSchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDACSchemeScalar_hh
