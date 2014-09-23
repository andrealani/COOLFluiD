#ifndef COOLFluiD_Numerics_FiniteElement_NumJacobStrategy_hh
#define COOLFluiD_Numerics_FiniteElement_NumJacobStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElement/ComputeJacobStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class BlockAccumulator; }

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy to compute the jacobian matrix of the system
 * via numerical diferenciation.
 *
 * @author Thomas Wuilbaut
 * @author Tiago Quintino
 * @author Pedro Maciel
 *
 */
class NumJacobStrategy : public ComputeJacobStrategy {
public:

  /**
   * Constructor.
   */
  NumJacobStrategy(const std::string& name);

  /**
   * Destructor.
   */
  ~NumJacobStrategy();

  /**
   * Add compute the term to add in the jacobian
   */
  void computeJacobianTerm();

  /**
   * Set up private data and data
   */
  void setup();

protected:

  /**
   * Set up private data and data
   */
  void cleanOtherResidual()
  {
    fill(_otherResidual.begin(),_otherResidual.end(),0.0);
  }

private:

  /// storage for the temporary node residuals
  RealVector _tempRes;

  /// storage for the temporary perturbed states
  std::vector<RealVector> _otherResidual;

}; // class NumJacobStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NumJacobStrategy_hh
