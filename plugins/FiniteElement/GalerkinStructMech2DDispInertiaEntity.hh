#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMech2DDispInertiaEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMech2DDispInertiaEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "InertiaEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "StructMech/StructMech2DInertiaDisp.hh"

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

class GalerkinStructMech2DDispInertiaEntity : public InertiaEntity {
public:

 /**
  * Default constructor
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  GalerkinStructMech2DDispInertiaEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMech2DDispInertiaEntity();

  /// Overloading of operator()
  RealMatrix& operator()();

  /// Setup the strategy
  virtual void setup();

private:

  /// Variable Set
  Common::SafePtr<Physics::StructMech::StructMech2DInertiaDisp> _structInertiaVarSet;

}; // end of class GalerkinStructMech2DDispInertiaEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMech2DDispInertiaEntity_hh

