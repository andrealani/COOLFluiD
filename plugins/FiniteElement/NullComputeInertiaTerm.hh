
#ifndef COOLFluiD_Numerics_FiniteElement_NullComputeInertiaTerm_hh
#define COOLFluiD_Numerics_FiniteElement_NullComputeInertiaTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "ComputeInertiaTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "InertiaEntity.hh"
#include "Common/SafePtr.hh"
#include "FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class Element;
  }

  namespace Numerics {

    namespace FiniteElement {


//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an object computing a Inertia term
 * with Finite Element method
 */
class NullComputeInertiaTerm :
  public ComputeInertiaTerm {


public:

  /// Constructor
  NullComputeInertiaTerm(const std::string& name);

  /// Default destructor
  ~NullComputeInertiaTerm();

  /// This ComputeTerm is not null
  bool isNull() const
  {
    return true;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NullComputeInertiaTerm";
  }

}; // end of class NullComputeInertiaTerm


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NullComputeInertiaTerm_hh

