#ifndef COOLFluiD_Fluctsplit_FluctSplitAdvectionDiffusion_hh
#define COOLFluiD_Fluctsplit_FluctSplitAdvectionDiffusion_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplit_API
/// @note build system defines FluctSplit_EXPORTS when compiling FluctSplit files
#ifdef FluctSplitAdvectionDiffusion_EXPORTS
#   define FluctSplitAdvectionDiffusion_API CF_EXPORT_API
#else
#   define FluctSplitAdvectionDiffusion_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitScalarDiffusion
class FluctSplitAdvectionDiffusionModule : public Environment::ModuleRegister<FluctSplitAdvectionDiffusionModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluctSplitAdvectionDiffusion";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() {  return "This module implementes Advection Diffusion model bindings to the Fluctuation Splitting space discretization method.";  }

}; // end FluctSplitNavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Fluctsplit_FluctSplitNavierStokes_hh
