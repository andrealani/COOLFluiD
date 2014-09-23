#ifndef COOLFluiD_Numerics_RemeshMeandr_hh
#define COOLFluiD_Numerics_RemeshMeandr_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement RemeshMeandr, a remesher based on
    /// system call to the Meandros mesh generation tool.
  namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module RemeshMeandr
class RemeshMeandrModule : public Environment::ModuleRegister<RemeshMeandrModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "RemeshMeandr";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements implements a remesher based on system call to the Meandros mesh generation tool.";
  }

}; // end RemeshMeandrModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandr

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_RemeshMeandr_hh

