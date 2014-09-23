#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinNeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinNeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "NeumannEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"

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

class GalerkinNeumannEntity : public NeumannEntity {

public:

  /// Default constructor without arguments
  GalerkinNeumannEntity(const std::string& name);

  /// Default destructor
  ~GalerkinNeumannEntity();

  /**
   * Setup the strategy
   */
  virtual void setup();

  /**
   * Overloading of operator()
   */
  virtual RealVector& operator()();

}; // end of class GalerkinNeumannEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinNeumannEntity_hh
