#ifndef COOLFluiD_FluxReconstructionNEQ_hh
#define COOLFluiD_FluxReconstructionNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Flux Reconstruction space discretization.
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionNEQ
 */
class FluxReconstructionNEQModule : public Environment::ModuleRegister<FluxReconstructionNEQModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionNEQ";
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

}; // end FluxReconstructionNEQModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_FluxReconstructionNEQ_hh

