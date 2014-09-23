#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinIndepSourceEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinIndepSourceEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "IndepSourceEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Independent
   * Source Term
   *
   * @author Tiago Quintino
   */
class GalerkinIndepSourceEntity : public IndepSourceEntity {
public:

  /// Default constructor without arguments
  GalerkinIndepSourceEntity(const std::string& name);

  /// Default destructor
  ~GalerkinIndepSourceEntity();

  /// Overloading of operator()
  RealVector& operator()();

}; // end of class GalerkinIndepSourceEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinIndepSourceEntity_hh

