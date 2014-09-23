#ifndef COOLFluiD_Numerics_FiniteElement_ComputeDiffusiveTerm_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeDiffusiveTerm_hh

#include "Common/NullableObject.hh"
#include "Framework/ComputeTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "DiffusiveEntity.hh"
#include "Common/SafePtr.hh"
#include "FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class DiffusiveVarSet;
    class Element;
  }

  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an object computing a diffusive term
 * with Finite Element method
 *
 * @author Tiago Quintino
 * @author Thomas Wuilbaut
 */
class ComputeDiffusiveTerm : public Framework::ComputeTerm<FiniteElementMethodData>{

public://types

  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ComputeDiffusiveTerm> PROVIDER;

public:

  /**
   * Constructor
   */
  ComputeDiffusiveTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeDiffusiveTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Setup the strategy
   */
  virtual void setup();


  /**
   * Is This ComputeTerm null
   */
  bool isNull() const
  {
    return false;
  }

  /**
   * Computes the Diffusive Term
   */
  void computeTerm(Framework::GeometricEntity* const cell, RealMatrix& result);

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ComputeDiffusiveTerm";
  }

protected:

  /// Diffusive Integrable Entity
  Common::SafePtr<DiffusiveEntity> _diffusiveEntity;

}; // end of class ComputeDiffusiveTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeDiffusiveTerm_hh
