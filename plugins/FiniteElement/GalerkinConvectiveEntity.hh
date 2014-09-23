#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinConvectiveEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinConvectiveEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConvectiveEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Convective Term
   */
class GalerkinConvectiveEntity : public ConvectiveEntity {

public:

  /// Default constructor without arguments
  GalerkinConvectiveEntity(const std::string& name);

  /// Default destructor
  ~GalerkinConvectiveEntity();

  /// Overloading of operator()
  virtual RealMatrix& operator()();

  /// Setup the strategy
  virtual void setup()
  {
    ConvectiveEntity::setup();
  }

}; // end of class GalerkinConvectiveEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinConvectiveEntity_hh
