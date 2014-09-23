#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinHeat2DPrimCoupledNeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinHeat2DPrimCoupledNeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "CoupledNeumannEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace Heat {
      class Heat2DDiffusivePrim;
    }
  }

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Diffusive Term
   *
   * @author Thomas Wuilbaut
   */

class GalerkinHeat2DPrimCoupledNeumannEntity : public CoupledNeumannEntity {

public:

  /// Default constructor without arguments
  GalerkinHeat2DPrimCoupledNeumannEntity(const std::string& name);

  /// Default destructor
  ~GalerkinHeat2DPrimCoupledNeumannEntity();

  /**
   * Setup the strategy
   */
  virtual void setup();

  /**
   * Overloading of operator()
   */
  virtual RealVector& operator()();

private:

   /// Variable Set
  Common::SafePtr<COOLFluiD::Physics::Heat::Heat2DDiffusivePrim> _heatDiffVarSet;

}; // end of class GalerkinHeat2DPrimCoupledNeumannEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinHeat2DPrimCoupledNeumannEntity_hh
