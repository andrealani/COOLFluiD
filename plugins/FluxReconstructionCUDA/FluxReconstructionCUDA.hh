#ifndef COOLFluiD_FluxReconstructionCUDA_hh
#define COOLFluiD_FluxReconstructionCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro FluxReconstruction_API
/// @note build system defines FluxReconstruction_EXPORTS
#ifdef FluxReconstructionCUDA_EXPORTS
#   define FluxReconstructionCUDA_API CF_EXPORT_API
#else
#   define FluxReconstructionCUDA_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module FluxReconstructionCUDA
class FluxReconstructionCUDAModule : public Environment::ModuleRegister<FluxReconstructionCUDAModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "FluxReconstructionCUDA"; }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements CUDA bindings for the Flux Reconstruction solver.";
  }
  
}; // end FluxReconstructionCUDAModule
      
//////////////////////////////////////////////////////////////////////////////

    }   // namespace FluxReconstruction

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionCUDA_hh
