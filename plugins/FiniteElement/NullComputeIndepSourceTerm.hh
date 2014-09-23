
#ifndef COOLFluiD_Numerics_FiniteElement_NullComputeIndepSourceTerm_hh
#define COOLFluiD_Numerics_FiniteElement_NullComputeIndepSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "ComputeIndepSourceTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "IndepSourceEntity.hh"
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
 * This class represents an object computing a IndepSource term
 * with Finite Element method
 */
class NullComputeIndepSourceTerm :
  public ComputeIndepSourceTerm {


public:

  /// Constructor
  NullComputeIndepSourceTerm(const std::string& name);

  /// Default destructor
  ~NullComputeIndepSourceTerm();

  /// This ComputeTerm is not null
  bool isNull() const
  {
    return true;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NullComputeIndepSourceTerm";
  }

}; // end of class NullComputeIndepSourceTerm


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NullComputeIndepSourceTerm_hh

