#ifndef COOLFluiD_Numerics_FiniteElement_ComputeLinearSourceTerm_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeLinearSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "Framework/ComputeTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "LinearSourceEntity.hh"
#include "Common/SafePtr.hh"
#include "FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class SourceVarSet;
    class Element;
  }

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an object computing a Source term
 * with Finite Element method
 *
 * @author Tiago Quintino
 * @author Thomas Wuilbaut
 *
 */
class ComputeLinearSourceTerm : public Framework::ComputeTerm<FiniteElementMethodData>{

public://types

  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ComputeLinearSourceTerm> PROVIDER;

public:

  /**
   * Constructor
   */
  ComputeLinearSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeLinearSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data to prepare the simulation
   */
  void setup();

  /**
   * Checks if the source term has linear part
   */
  bool hasLinearCoef()
  {
    return getMethodData().getSourceVar()->hasLinearCoef();
  }

  /**
   * This ComputeTerm is null if the SourceVarSet is null
   */
  bool isNull() const
  {
    return false;
  }

  /**
   * Computes the Source Term
   */
  void computeTerm(Framework::GeometricEntity* const cell, RealMatrix& result);

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ComputeLinearSourceTerm";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected:

  /// Source Integrable Entity
  Common::SafePtr<LinearSourceEntity> _linearSourceEntity;

}; // end of class ComputeLinearSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeLinearSourceTerm_hh
