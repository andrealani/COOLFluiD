#ifndef COOLFluiD_FluxReconstructionMultiFluidMHD_hh
#define COOLFluiD_FluxReconstructionMultiFluidMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Flux Reconstruction space discretization.
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionMultiFluidMHD
 */
class FluxReconstructionMultiFluidMHDModule : public Environment::ModuleRegister<FluxReconstructionMultiFluidMHDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionMultiFluidMHD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Flux Reconstruction space discretization for Navier-Stokes.";
  }

}; // end FluxReconstructionMultiFluidMHDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_FluxReconstructionMultiFluidMHD_hh

