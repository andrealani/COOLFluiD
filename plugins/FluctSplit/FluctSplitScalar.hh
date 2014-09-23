#ifndef COOLFluiD_Numerics_Fluctsplit_FluctSplitScalar_hh
#define COOLFluiD_Numerics_Fluctsplit_FluctSplitScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluctSplitScalar_API
/// @note build system defines FluctSplitScalar_EXPORTS
#ifdef FluctSplitScalar_EXPORTS
#   define FluctSplitScalar_API CF_EXPORT_API
#else
#   define FluctSplitScalar_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement Fluctuation Splitting space method.
    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluctSplitScalar
class FluctSplitScalarModule : public Environment::ModuleRegister<FluctSplitScalarModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "FluctSplitScalar";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes the scalar schemes of the Fluctuation Splitting space discretization method.";
  }

}; // end FluctSplitScalarModule

//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Fluctsplit_FluctSplitScalar_hh
