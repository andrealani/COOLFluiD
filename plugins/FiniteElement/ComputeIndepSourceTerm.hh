#ifndef COOLFluiD_Numerics_FiniteElement_ComputeIndepSourceTerm_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeIndepSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "Framework/ComputeTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "IndepSourceEntity.hh"
#include "Common/SafePtr.hh"
#include "Framework/VectorialFunction.hh"
#include "FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
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
class ComputeIndepSourceTerm : public Framework::ComputeTerm<FiniteElementMethodData>{

public://types

  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ComputeIndepSourceTerm> PROVIDER;


public:

  /**
   * Constructor
   */
  ComputeIndepSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeIndepSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * This ComputeTerm is null if the SourceVarSet is null
   */
  bool isNull() const
  {
    return false;
  }

  /**
   * Checks if the source term has linear part
   */
  bool hasIndepCoef()
  {
    return getMethodData().getSourceVar()->hasIndepCoef();
  }

  /**
   * Set up private data to prepare the simulation
   */
  void setup();

  /**
   * Computes the Source Term
   */
  void computeTerm(Framework::GeometricEntity* const cell, RealVector& result);

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ComputeIndepSourceTerm";
  }
 
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected:

  /// Source Integrable Entity
  Common::SafePtr<IndepSourceEntity> _indepSourceEntity;

}; // end of class ComputeIndepSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeIndepSourceTerm_hh
