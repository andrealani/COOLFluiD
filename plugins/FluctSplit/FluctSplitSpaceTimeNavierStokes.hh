#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitSpaceTimeNavierStokes_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitSpaceTimeNavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitSpaceTimeNavierStokes_API
/// @note build system defines FluctSplitSpaceTimeNavierStokes_EXPORTS
#ifdef FluctSplitSpaceTimeNavierStokes_EXPORTS
#   define FluctSplitSpaceTimeNavierStokes_API CF_EXPORT_API
#else
#   define FluctSplitSpaceTimeNavierStokes_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitSpaceTimeNavierStokes
class FluctSplitSpaceTimeNavierStokesModule : public Environment::ModuleRegister<FluctSplitSpaceTimeNavierStokesModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitSpaceTimeNavierStokes";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements the space time schemes of the Fluctuation Splitting space discretization method for Navier-Stokes.";
  }

}; // end FluctSplitSpaceTimeNavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitSpaceTimeNavierStokes_hh
