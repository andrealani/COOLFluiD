#ifndef COOLFluiD_Numerics_FiniteElement_IndepSourceEntity_hh
#define COOLFluiD_Numerics_FiniteElement_IndepSourceEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEMIntegrableEntity.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/SafePtr.hh"
#include "Framework/SourceVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity fo a diffusive term
   *
   * @author Tiago Quintino
   *
   * @see IntegrableEntity
   *
   */
class IndepSourceEntity : public FiniteElement::FEMIntegrableEntity {
public:

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,IndepSourceEntity> PROVIDER;

 /**
  * Default constructor without arguments
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  IndepSourceEntity(const std::string& name);

 /**
  * Default destructor
  */
  ~IndepSourceEntity();

  /**
   * Set the Source variable set
   */
  void setup()
  {
    FEMIntegrableEntity::setup();

    _sourceVarSet = getMethodData().getSourceVar();
    _nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    _coefs.resize(_nbEqs);
    _result.resize(_nbEqs);
  }

  /**
   * Overloading of operator()
   */
  virtual RealVector& operator()() = 0;

  /**
   * Set the size of the problem
   */
  CFuint size()
  {
    return _result.size();
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "IndepSourceEntity";
  }

protected:

  /// Source Variable Set
  Common::SafePtr<Framework::SourceVarSet> _sourceVarSet;

  /// Diffusive Coefficients Matrix Storage
  RealVector _coefs;

  /// vector storing the temporary computation
  RealVector _result;

  ///  nbEqs
  CFuint _nbEqs;

}; // end of class IndepSourceEntity

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteElement

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_IndepSourceEntity_hh
