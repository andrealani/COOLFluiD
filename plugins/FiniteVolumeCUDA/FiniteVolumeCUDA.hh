#ifndef COOLFluiD_Numerics_FiniteVolume_FiniteVolumeCUDA_hh
#define COOLFluiD_Numerics_FiniteVolume_FiniteVolumeCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FiniteVolume_API
/// @note build system defines FiniteVolume_EXPORTS
#ifdef FiniteVolumeCUDA_EXPORTS
#   define FiniteVolumeCUDA_API CF_EXPORT_API
#else
#   define FiniteVolumeCUDA_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FiniteVolumeCUDA
class FiniteVolumeCUDAModule : public Environment::ModuleRegister<FiniteVolumeCUDAModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FiniteVolumeCUDA"; }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements CUDA bindings for the Finite Volume solver.";
  }
  
}; // end FiniteVolumeCUDAModule
      
//////////////////////////////////////////////////////////////////////////////

    }   // namespace FiniteVolume

  }   // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FiniteVolumeCUDA_hh
