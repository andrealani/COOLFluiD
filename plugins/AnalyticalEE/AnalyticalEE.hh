#ifndef COOLFluiD_Numerics_AnalyticalEE_hh
#define COOLFluiD_Numerics_AnalyticalEE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement an Error estimator which compares the computed
  /// solution with the a analytical expression at each State
  namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module AnalyticalEE
/// @author Tiago Quintino
class AnalyticalEEModule : public Environment::ModuleRegister<AnalyticalEEModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "AnalyticalEE";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements an Error estimator which compares the computed"
           " solution qith the a analytical expression at each State";
  }

}; // end AnalyticalEEModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_AnalyticalEE_hh

