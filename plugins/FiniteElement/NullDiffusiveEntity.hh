#ifndef COOLFluiD_Numerics_FiniteElement_NullDiffusiveEntity_hh
#define COOLFluiD_Numerics_FiniteElement_NullDiffusiveEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiffusiveEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Null IntegrableEntity for a Diffusive Term
   *
   * @author Tiago Quintino
   */
class NullDiffusiveEntity : public DiffusiveEntity {

public:

  /// Default constructor without arguments
  NullDiffusiveEntity(const std::string& name);

  /// Default destructor
  virtual ~NullDiffusiveEntity();

  /**
   * Overloading of operator()
   */
  virtual RealMatrix& operator()();

}; // end of class NullDiffusiveEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NullDiffusiveEntity_hh
