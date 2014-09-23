#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitSystem_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitSystem_API
/// @note build system defines FluctSplitSystem_EXPORTS
#ifdef FluctSplitSystem_EXPORTS
#   define FluctSplitSystem_API CF_EXPORT_API
#else
#   define FluctSplitSystem_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitSystem
class FluctSplitSystemModule : public Environment::ModuleRegister<FluctSplitSystemModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitSystem";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes the system schemes of the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitSystemModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitSystem_hh
