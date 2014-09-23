#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinHeat3DPrimDiffEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinHeat3DPrimDiffEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiffusiveEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace Heat {
      class Heat3DDiffusivePrim;
    }
  }

  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Diffusive Term
   *
   * @author Tiago Quintino
   */
class GalerkinHeat3DPrimDiffEntity : public DiffusiveEntity {

public:

  /// Default constructor without arguments
  GalerkinHeat3DPrimDiffEntity(const std::string& name);

  /// Default destructor
  virtual ~GalerkinHeat3DPrimDiffEntity();

  /**
   * Setup the strategy
   */
  virtual void setup();

  /**
   * Overloading of operator()
   */
  virtual RealMatrix& operator()();

private:

  /// Variable Set
  Common::SafePtr<COOLFluiD::Physics::Heat::Heat3DDiffusivePrim> _heatDiffVarSet;

}; // end of class GalerkinHeat3DPrimDiffEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinHeat3DPrimDiffEntity_hh
