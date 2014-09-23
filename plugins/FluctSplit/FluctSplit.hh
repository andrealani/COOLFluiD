#ifndef COOLFluiD_FluctSplit_hh
#define COOLFluiD_FluctSplit_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplit_API
/// @note build system defines FluctSplit_EXPORTS when compiling FluctSplit files
#ifdef FluctSplit_EXPORTS
#   define FluctSplit_API CF_EXPORT_API
#else
#   define FluctSplit_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplit
class FluctSplit_API FluctSplitModule : public Environment::ModuleRegister<FluctSplitModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluctSplit";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() {  return "This module implementes the Fluctuation Splitting space discretization method."; }

}; // end FluctSplitModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_hh
