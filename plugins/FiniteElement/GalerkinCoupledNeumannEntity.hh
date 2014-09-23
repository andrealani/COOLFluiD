#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinCoupledNeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinCoupledNeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "CoupledNeumannEntity.hh"
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

class GalerkinCoupledNeumannEntity : public CoupledNeumannEntity {

public:

  /// Default constructor without arguments
  GalerkinCoupledNeumannEntity(const std::string& name);

  /// Default destructor
  ~GalerkinCoupledNeumannEntity();

  /**
   * Setup the strategy
   */
  virtual void setup();

  /**
   * Overloading of operator()
   */
  virtual RealVector& operator()();

}; // end of class GalerkinCoupledNeumannEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinCoupledNeumannEntity_hh
