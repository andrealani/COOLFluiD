#ifndef COOLFluiD_Numerics_FluctSplit_STM_LDASchemeScalarHO_hh
#define COOLFluiD_Numerics_FluctSplit_STM_LDASchemeScalarHO_hh

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
class STM_LDASchemeScalarHO : public STM_SplitterScalar {
public:

  /**
   * Default constructor.
   */
  STM_LDASchemeScalarHO(const std::string& name);

  /**
   * Default destructor
   */
  ~STM_LDASchemeScalarHO();

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
  void distributePast(std::vector<RealVector>& residual);

  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealVector _sumKplus;

}; // end of class STM_LDASchemeScalarHO

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeScalar_hh
