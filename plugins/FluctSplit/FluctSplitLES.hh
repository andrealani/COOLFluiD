#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitLES_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitLES_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitNavierStokes_API
/// @note build system defines FluctSplitNavierStokes_EXPORTS
#ifdef FluctSplitLES_EXPORTS
#   define FluctSplitLES_API CF_EXPORT_API
#else
#   define FluctSplitLES_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitLES
class FluctSplitLESModule : public Environment::ModuleRegister<FluctSplitLESModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitLES";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes LES model bindings to the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitLESModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitNavierStokes_hh
