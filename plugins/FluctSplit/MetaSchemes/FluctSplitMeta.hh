#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitMeta_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitMeta_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitMeta_API
/// @note build system defines FluctSplitMeta_EXPORTS when compiling FluctSplit files
#ifdef FluctSplitMeta_EXPORTS
#   define FluctSplitMeta_API CF_EXPORT_API
#else
#   define FluctSplitMeta_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitMeta
class FluctSplitMetaModule : public Environment::ModuleRegister<FluctSplitMetaModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluctSplitMeta";   }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes meta schemes for the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitMetaModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitMeta_hh
