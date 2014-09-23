#ifndef COOLFluiD_Numerics_FiniteElement_ConvectiveEntity_hh
#define COOLFluiD_Numerics_FiniteElement_ConvectiveEntity_hh

#include "Framework/GeometricEntity.hh"
#include "FEMIntegrableEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/SafePtr.hh"
#include "Framework/ConvectiveVarSet.hh"

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {


  /**
   * This class represents a IntegrableEntity for a convective term
   *
   * @see IntegrableEntity
   */
class ConvectiveEntity : public FEMIntegrableEntity {
public: // methods

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ConvectiveEntity> PROVIDER;

  /// Default constructor
  ConvectiveEntity(const std::string& name);

  /// Default destructor
  ~ConvectiveEntity();

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

    _convVarSet = getMethodData().getUpdateVar();
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
    return "ConvectiveEntity";
  }

protected: // data

  /// Convective Variable Set
  Common::SafePtr<Framework::ConvectiveVarSet> _convVarSet;

  /// Convective Coefficients Matrix Storage
  RealMatrix _coefs;

  /// matrix storing the temporary computation
  RealMatrix _result;

  /// number of equations
  CFuint _nbEqs;

}; // end of class ConvectiveEntity


    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

#endif // COOLFluiD_Numerics_FiniteElement_ConvectiveEntity_hh

