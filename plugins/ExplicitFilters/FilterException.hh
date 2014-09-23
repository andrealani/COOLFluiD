#ifndef COOLFluiD_Numerics_ExplicitFilters_FilterException_hh
#define COOLFluiD_Numerics_ExplicitFilters_FilterException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace ExplicitFilters {
  
//////////////////////////////////////////////////////////////////////////////

/// Exceptions thrown within the ExplicitFilters plugin
class Framework_API FilterException : public Common::Exception {
public:

  /// Constructor
  FilterException(const Common::CodeLocation& where, const std::string& what, const std::string& className)
    : Common::Exception(where, what, className) {}
     
  FilterException(const Common::CodeLocation& where, const std::string& what)
    : Common::Exception(where, what, "ExplicitFilter Exception") {}
    
  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  FilterException(const FilterException& e) throw() : Common::Exception(e) {}
    
}; // end of class FilterException

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_FilterException_hh
