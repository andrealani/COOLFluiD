#ifndef COOLFluiD_FluxReconstructionPoisson_hh
#define COOLFluiD_FluxReconstructionPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Flux Reconstruction space discretization.
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionPoisson
 */
class FluxReconstructionPoissonModule : public Environment::ModuleRegister<FluxReconstructionPoissonModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionPoisson";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Flux Reconstruction space discretization for Poisson.";
  }

}; // end FluxReconstructionPoissonModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_FluxReconstructionPoisson_hh

