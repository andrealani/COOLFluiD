#ifndef COOLFluiD_Numerics_FiniteElement_NullComputeConvectiveTerm_hh
#define COOLFluiD_Numerics_FiniteElement_NullComputeConvectiveTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "ComputeConvectiveTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "ConvectiveEntity.hh"
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
 * This class represents an object computing a convective term
 * with Finite Element method
 */
class NullComputeConvectiveTerm :
  public ComputeConvectiveTerm {


public:

  /// Constructor
  NullComputeConvectiveTerm(const std::string& name);

  /// Default destructor
  ~NullComputeConvectiveTerm();

  /// This ComputeTerm is not null
  bool isNull() const
  {
    return true;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NullComputeConvectiveTerm";
  }

}; // end of class NullComputeConvectiveTerm


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NullComputeConvectiveTerm_hh

