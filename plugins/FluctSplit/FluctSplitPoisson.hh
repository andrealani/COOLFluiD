#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitPoisson_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitPoisson_API
/// @note build system defines FluctSplitPoisson_EXPORTS
#ifdef FluctSplitPoisson_EXPORTS
#   define FluctSplitPoisson_API CF_EXPORT_API
#else
#   define FluctSplitPoisson_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitPoisson
class FluctSplitPoissonModule : public Environment::ModuleRegister<FluctSplitPoissonModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitPoisson";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes Poisson model bindings to the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitPoissonModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitPoisson_hh
