#ifndef COOLFluiD_Numerics_FiniteElement_LinearSourceEntity_hh
#define COOLFluiD_Numerics_FiniteElement_LinearSourceEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEMIntegrableEntity.hh"
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
class LinearSourceEntity : public FiniteElement::FEMIntegrableEntity {
public:

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,LinearSourceEntity> PROVIDER;

 /**
  * Default constructor without arguments
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  LinearSourceEntity(const std::string& name);

 /**
  * Default destructor
  */
  ~LinearSourceEntity();

  /**
   * Set the Source variable set
   */
  void setup()
  {
    FEMIntegrableEntity::setup();

    _sourceVarSet = getMethodData().getSourceVar();
    _nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    _coefs.resize(_nbEqs,_nbEqs);
    _result.resize(_nbEqs,_nbEqs);
  }

 /**
  * @return the rows of the result
  */
  CFuint getNbRows() const
  {
   return _result.nbRows();
  }

 /**
  * @return the number of columns of the result
  */
  CFuint getNbColumns() const
  {
    return _result.nbCols();
  }

  /**
   * Overloading of operator()
   */
  virtual RealMatrix& operator()() = 0;

  /**
   * Get the size of the result
   */
  CFuint size() { return _result.size(); }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "LinearSourceEntity";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

protected:

  /// Source Variable Set
  Common::SafePtr<Framework::SourceVarSet> _sourceVarSet;

  /// Diffusive Coefficients Matrix Storage
  RealMatrix _coefs;

  /// matrix storing the temporary computation
  RealMatrix _result;

  ///  state i
  CFuint _nbEqs;

}; // end of class LinearSourceEntity

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteElement

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_LinearSourceEntity_hh
