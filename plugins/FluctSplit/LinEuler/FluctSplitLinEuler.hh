#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitLinEuler_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitLinEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitLinEuler_API
/// @note build system defines FluctSplitLinEuler_EXPORTS
#ifdef FluctSplitLinEuler_EXPORTS
#   define FluctSplitLinEuler_API CF_EXPORT_API
#else
#   define FluctSplitLinEuler_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitLinEuler
class FluctSplitLinEuler_API FluctSplitLinEulerModule : public Environment::ModuleRegister<FluctSplitLinEulerModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluctSplitLinEuler";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() { return "This module implementes Linearized Euler model bindings to the Fluctuation Splitting space discretization method."; }

}; // end FluctSplitLinEulerModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitLinEuler_hh
