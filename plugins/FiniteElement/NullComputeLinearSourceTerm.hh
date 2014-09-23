
#ifndef COOLFluiD_Numerics_FiniteElement_NullComputeLinearSourceTerm_hh
#define COOLFluiD_Numerics_FiniteElement_NullComputeLinearSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "ComputeLinearSourceTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "LinearSourceEntity.hh"
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
 * This class represents an object computing a LinearSource term
 * with Finite Element method
 */
class NullComputeLinearSourceTerm :
  public ComputeLinearSourceTerm {


public:

  /// Constructor
  NullComputeLinearSourceTerm(const std::string& name);

  /// Default destructor
  ~NullComputeLinearSourceTerm();

  /// This ComputeTerm is not null
  bool isNull() const
  {
    return true;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NullComputeLinearSourceTerm";
  }

}; // end of class NullComputeLinearSourceTerm


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NullComputeLinearSourceTerm_hh

