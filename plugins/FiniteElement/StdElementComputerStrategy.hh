#ifndef COOLFluiD_Numerics_FiniteElement_StdElementComputerStrategy_hh
#define COOLFluiD_Numerics_FiniteElement_StdElementComputerStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeTerm.hh"

#include "FiniteElement/ComputeResidualStrategy.hh"
#include "FiniteElement/FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a strategy to compute the residual
 *
 * @author Thomas Wuilbaut
 * @author Tiago Quintino
 * @author Pedro Maciel
 *
 */
class StdElementComputerStrategy : public ComputeResidualStrategy {
public:

  /**
   * Constructor.
   */
  StdElementComputerStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdElementComputerStrategy();

  /**
   * Computes the element matrix
   */
  void computeElemMatrix();

  /**
   * Computes the element vector
   */
  void computeElemVector();

  /**
   * Set up private data and data
   */
  void setup();

private: // data

  /// Temporary storage of the integration result
  RealMatrix _integResultMat;

  /// Temporary storage of the integration result
  RealVector _integResultVec;

  /// list of term computers that have a matrix has result
  std::vector<Common::SafePtr<Framework::ComputeTerm<FiniteElementMethodData> > > _matrixComputeTerms;

  /// list of term computers that have a vector has result
  std::vector<Common::SafePtr<Framework::ComputeTerm<FiniteElementMethodData> > > _vectorComputeTerms;

}; // class StdElementComputerStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_StdElementComputerStrategy_hh
