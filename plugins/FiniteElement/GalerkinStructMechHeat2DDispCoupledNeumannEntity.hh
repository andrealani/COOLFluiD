#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DDispCoupledNeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DDispCoupledNeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "CoupledNeumannEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"
#include "StructMechHeat/StructMechHeat2DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Diffusive Term
   *
   * @author Thomas Wuilbaut
   */

class GalerkinStructMechHeat2DDispCoupledNeumannEntity : public CoupledNeumannEntity {

public:

  /// Default constructor without arguments
  GalerkinStructMechHeat2DDispCoupledNeumannEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMechHeat2DDispCoupledNeumannEntity();

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
  Common::SafePtr<Physics::StructMechHeat::StructMechHeat2DDiffusiveVarSet> _structDiffVarSet;

}; // end of class GalerkinStructMechHeat2DDispCoupledNeumannEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DDispCoupledNeumannEntity_hh
