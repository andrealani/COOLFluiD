#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitHONavierStokes_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitHONavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplit_API
/// @note build system defines FluctSplit_EXPORTS
#ifdef FluctSplitHONavierStokes_EXPORTS
#   define FluctSplitHONavierStokes_API CF_EXPORT_API
#else
#   define FluctSplitHONavierStokes_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitHONavierStokes
class FluctSplitHONavierStokesModule : public Environment::ModuleRegister<FluctSplitHONavierStokesModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluctSplitHONavierStokes"; }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes NavierStokes model bindings to the High Order Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitHONavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitHONavierStokes_hh
