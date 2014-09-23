#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitMHD_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitMHD_API
/// @note build system defines FluctSplitMHD_EXPORTS
#ifdef FluctSplitMHD_EXPORTS
#   define FluctSplitMHD_API CF_EXPORT_API
#else
#   define FluctSplitMHD_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitMHD
class FluctSplitMHDModule : public Environment::ModuleRegister<FluctSplitMHDModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitMHD";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes MHD model bindings to the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitMHDModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitMHD_hh
