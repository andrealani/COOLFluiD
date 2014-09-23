#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitKOmega_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitKOmega_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitKOmega_API
/// @note build system defines FluctSplitKOmega_EXPORTS
#ifdef FluctSplitKOmega_EXPORTS
#   define FluctSplitKOmega_API CF_EXPORT_API
#else
#   define FluctSplitKOmega_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitKOmega
class FluctSplitKOmega_API FluctSplitKOmegaModule : public Environment::ModuleRegister<FluctSplitKOmegaModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluctSplitKOmega";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() { return "This module implementes k-omega RANS model bindings to the Fluctuation Splitting space discretization method."; }

}; // end FluctSplitKOmegaModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitKOmega_hh
