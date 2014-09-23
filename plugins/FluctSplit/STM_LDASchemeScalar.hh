#ifndef COOLFluiD_Numerics_FluctSplit_STM_LDASchemeScalar_hh
#define COOLFluiD_Numerics_FluctSplit_STM_LDASchemeScalar_hh

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
class STM_LDASchemeScalar : public STM_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STM_LDASchemeScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_LDASchemeScalar();

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

}; // end of class STM_LDASchemeScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeScalar_hh
