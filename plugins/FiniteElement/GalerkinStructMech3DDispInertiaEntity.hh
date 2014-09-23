#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMech3DDispInertiaEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMech3DDispInertiaEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "InertiaEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "StructMech/StructMech3DInertiaDisp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Inertia Term
   *
   * @author Thomas Wuilbaut
   */

class GalerkinStructMech3DDispInertiaEntity : public InertiaEntity {
public:

 /**
  * Default constructor
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  GalerkinStructMech3DDispInertiaEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMech3DDispInertiaEntity();

  /// Overloading of operator()
  RealMatrix& operator()();

  /// Setup the strategy
  virtual void setup();

private:

  /// Variable Set
  Common::SafePtr<Physics::StructMech::StructMech3DInertiaDisp> _structInertiaVarSet;

}; // end of class GalerkinStructMech3DDispInertiaEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMech3DDispInertiaEntity_hh

