#ifndef COOLFluiD_FluxReconstructionMHD_hh
#define COOLFluiD_FluxReconstructionMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Flux Reconstruction space discretization.
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionMHD
 */
class FluxReconstructionMHDModule : public Environment::ModuleRegister<FluxReconstructionMHDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionMHD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Flux Reconstruction space discretization for MHD.";
  }

}; // end FluxReconstructionMHDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_FluxReconstructionMHD_hh

