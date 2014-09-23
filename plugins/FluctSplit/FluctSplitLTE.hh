#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitLTE_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitLTE_API
/// @note build system defines FluctSplitLTE_EXPORTS
#ifdef FluctSplitLTE_EXPORTS
#   define FluctSplitLTE_API CF_EXPORT_API
#else
#   define FluctSplitLTE_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitLTE
class FluctSplitLTEModule : public Environment::ModuleRegister<FluctSplitLTEModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitLTE";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes NavierStokes model bindings to the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitLTEModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitLTE_hh
