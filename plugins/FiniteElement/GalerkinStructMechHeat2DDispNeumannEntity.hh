#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DDispNeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DDispNeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "NeumannEntity.hh"
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

class GalerkinStructMechHeat2DDispNeumannEntity : public NeumannEntity {

public:

  /// Default constructor without arguments
  GalerkinStructMechHeat2DDispNeumannEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMechHeat2DDispNeumannEntity();

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

}; // end of class GalerkinStructMechHeat2DDispNeumannEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMechHeat2DDispNeumannEntity_hh
