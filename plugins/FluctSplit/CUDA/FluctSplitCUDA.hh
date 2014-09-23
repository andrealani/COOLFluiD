#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitCUDA_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplit_API
/// @note build system defines FluctSplit_EXPORTS
#ifdef FluctSplitCUDA_EXPORTS
#   define FluctSplitCUDA_API CF_EXPORT_API
#else
#   define FluctSplitCUDA_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitCUDA
class FluctSplitCUDAModule : public Environment::ModuleRegister<FluctSplitCUDAModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluctSplitCUDA"; }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes NavierStokes model bindings to the High Order Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitCUDAModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitCUDA_hh
