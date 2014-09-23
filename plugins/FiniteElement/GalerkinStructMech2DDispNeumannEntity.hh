#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinStructMech2DDispNeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinStructMech2DDispNeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "NeumannEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"
#include "StructMech/StructMech2DDiffusiveVarSet.hh"

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

class GalerkinStructMech2DDispNeumannEntity : public NeumannEntity {

public:

  /// Default constructor without arguments
  GalerkinStructMech2DDispNeumannEntity(const std::string& name);

  /// Default destructor
  ~GalerkinStructMech2DDispNeumannEntity();

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
  Common::SafePtr<Physics::StructMech::StructMech2DDiffusiveVarSet> _structDiffVarSet;

}; // end of class GalerkinStructMech2DDispNeumannEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinStructMech2DDispNeumannEntity_hh
