#ifndef COOLFluiD_Numerics_FiniteElement_DiffusiveEntity_hh
#define COOLFluiD_Numerics_FiniteElement_DiffusiveEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "Common/SafePtr.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "FiniteElement/FEMIntegrableEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a IntegrableEntity for a diffusive term
 *
 * @author Tiago Quintino
 *
 * @see IntegrableEntity
 */
class DiffusiveEntity : public FiniteElement::FEMIntegrableEntity {
public:

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,DiffusiveEntity> PROVIDER;

  /// Default constructor
  DiffusiveEntity(const std::string& name);

  /// Default destructor
  virtual ~DiffusiveEntity();

  /// @return the rows of the result
  CFuint getNbRows() const
  {
    return _result.nbRows();
  }

  /// @return the number of columns of the result
  CFuint getNbColumns() const
  {
    return _result.nbCols();
  }

  /**
   * Overloading of operator()
   * @pre must have first setVarSet()
   */
  virtual RealMatrix& operator()() = 0;

  /// Setup the strategy
  virtual void setup()
  {
    FEMIntegrableEntity::setup();

    _diffVarSet = getMethodData().getDiffusiveVar();
    _nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    _result.resize(_nbEqs,_nbEqs);
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
    return "DiffusiveEntity";
  }

protected:

  /// Diffusive Variable Set
  Common::SafePtr<Framework::DiffusiveVarSet> _diffVarSet;

  /// Diffusive Coefficients Matrix Storage
  RealMatrix _coefs;

  /// matrix storing the temporary computation
  RealMatrix _result;

  ///  number of equations
  CFuint _nbEqs;

}; // end of class DiffusiveEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_DiffusiveEntity_hh
