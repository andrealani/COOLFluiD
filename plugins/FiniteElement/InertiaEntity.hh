#ifndef COOLFluiD_Numerics_FiniteElement_InertiaEntity_hh
#define COOLFluiD_Numerics_FiniteElement_InertiaEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEMIntegrableEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/SafePtr.hh"
#include "Framework/InertiaVarSet.hh"

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
class InertiaEntity : public FiniteElement::FEMIntegrableEntity {
public: // methods

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,InertiaEntity> PROVIDER;

 /**
  * Default constructor without arguments
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  InertiaEntity(const std::string& name);

 /**
  * Default destructor
  */
  ~InertiaEntity();

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
   * Set the size of the problem
   */
  virtual void setup()
  {
    FEMIntegrableEntity::setup();

    _inertiaVarSet = getMethodData().getInertiaVar();
    _nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    _result.resize(_nbEqs, _nbEqs);
  }

  /**
   * Get the size of the result
   */
  CFuint size() { return _result.size(); }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "InertiaEntity";
  }

protected: // data

  /// Inertia Variable Set
  Common::SafePtr<Framework::InertiaVarSet> _inertiaVarSet;

  /// matrix storing the temporary computation
  RealMatrix _result;

  ///  state i
  CFuint _nbEqs;

}; // end of class InertiaEntity

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteElement

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_InertiaEntity_hh
