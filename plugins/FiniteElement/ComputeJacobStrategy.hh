#ifndef COOLFluiD_Numerics_FiniteElement_ComputeJacobStrategy_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeJacobStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntity.hh"
#include "Framework/BaseMethodStrategyProvider.hh"

#include "FiniteElement/FiniteElementMethodData.hh"
#include "FiniteElement/FElemTypeData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class BlockAccumulator; }

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a finite element strategy to compute a jacobian element matrix
/// @author Thomas Wuilbaut
/// @author Tiago Quintino
/// @author Pedro Maciel
class ComputeJacobStrategy : public FiniteElementMethodStrategy {
public: // types

  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ComputeJacobStrategy> PROVIDER;

public: // methods

  /// Constructor.
  ComputeJacobStrategy(const std::string& name);

  /// Destructor.
  virtual ~ComputeJacobStrategy();

  /// Add compute the term to add in the jacobian
  virtual void computeJacobianTerm() = 0;

  /// Gets the Class name
  static std::string getClassName() { return "ComputeJacobStrategy"; }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /// Set up private data and data
  virtual void setup();

  /// Set up private data and data
  virtual void unsetup();

}; // class ComputeJacobStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeJacobStrategy_hh
