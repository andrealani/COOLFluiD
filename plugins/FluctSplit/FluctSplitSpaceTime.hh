#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitSpaceTime_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitSpaceTime_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitSpaceTime_API
/// @note build system defines FluctSplitSpaceTime_EXPORTS
#ifdef FluctSplitSpaceTime_EXPORTS
#   define FluctSplitSpaceTime_API CF_EXPORT_API
#else
#   define FluctSplitSpaceTime_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitSpaceTime
class FluctSplitSpaceTimeModule : public Environment::ModuleRegister<FluctSplitSpaceTimeModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitSpaceTime";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes the space time schemes of the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitSpaceTimeModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitSpaceTime_hh
