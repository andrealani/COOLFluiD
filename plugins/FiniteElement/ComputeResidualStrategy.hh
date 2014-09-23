#ifndef COOLFluiD_Numerics_FiniteElement_ComputeResidualStrategy_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeResidualStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/BlockAccumulator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the element residual of the form:
 *
 *  R(U) = K*U -f
 *
 * @author Thomas Wuilbaut
 * @author Tiago Quintino
 * @author Pedro Maciel
 *
 */
class ComputeResidualStrategy : public FiniteElementMethodStrategy {
public: // types

  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ComputeResidualStrategy> PROVIDER;

public: // methods

  /**
   * Constructor.
   */
  ComputeResidualStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ComputeResidualStrategy();

  /**
   * Computes the residual
   */
  virtual void computeElementResidual(std::vector<RealVector>& residual);

  /**
   * Computes the element matrix
   * @param elemMat square matrix sized nbStates*nbEquations which will be overwritten. No need to clear before.
   */
  virtual void computeElemMatrix() = 0;

  /**
   * Computes the element vector
   * @param elemVec vector sized nbStates*nbEquations which will be overwritten. No need to clear before.
   */
  virtual void computeElemVector() = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ComputeResidualStrategy";
  }

  /**
   * Set up private data and data
   */
  virtual void setup();

protected: // data

}; // class ComputeResidualStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeResidualStrategy_hh
