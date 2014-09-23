#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinLinearSourceEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinLinearSourceEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "LinearSourceEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Linear Source Term
   *
   * @author Tiago Quintino
   */
class GalerkinLinearSourceEntity : public LinearSourceEntity {
public:

  /// Default constructor without arguments
  GalerkinLinearSourceEntity(const std::string& name);

  /// Default destructor
  ~GalerkinLinearSourceEntity();

  /// Overloading of operator()
  RealMatrix& operator()();

}; // end of class GalerkinLinearSourceEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinLinearSourceEntity_hh

