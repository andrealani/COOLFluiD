#ifndef COOLFluiD_Numerics_FiniteElement_NullInertiaEntity_hh
#define COOLFluiD_Numerics_FiniteElement_NullInertiaEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "InertiaEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Null IntegrableEntity for a Inertia Term
   *
   * @author Thomas Wuilbaut
   */

class NullInertiaEntity : public InertiaEntity {
public:

 /**
  * Default constructor
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  NullInertiaEntity(const std::string& name);

  /// Default destructor
  ~NullInertiaEntity();

  /// Overloading of operator()
  RealMatrix& operator()();

}; // end of class NullInertiaEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NullInertiaEntity_hh

