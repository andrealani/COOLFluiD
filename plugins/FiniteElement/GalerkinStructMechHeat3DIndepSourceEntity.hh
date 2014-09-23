#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat3DIndepSourceEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat3DIndepSourceEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "IndepSourceEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealVector.hh"
#include "StructMechHeat/StructMechHeat3DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Independent
   * Source Term
   *
   * @author Thomas Wuilbaut
   */
class GalerkinStructMechHeat3DIndepSourceEntity : public IndepSourceEntity {
public:

  /// Default constructor without arguments
  GalerkinStructMechHeat3DIndepSourceEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMechHeat3DIndepSourceEntity();

  /// Setup
  void setup();

  /// Overloading of operator()
  RealVector& operator()();

protected:

  /// Variable Set
  Common::SafePtr<Physics::StructMechHeat::StructMechHeat3DDiffusiveVarSet> _structDiffVarSet;

}; // end of class GalerkinStructMechHeat3DIndepSourceEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat3DIndepSourceEntity_hh

