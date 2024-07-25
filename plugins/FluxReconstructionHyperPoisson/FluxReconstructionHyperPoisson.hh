#ifndef COOLFluiD_FluxReconstructionHyperPoisson_hh
#define COOLFluiD_FluxReconstructionHyperPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Flux Reconstruction space discretization.
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionHyperPoisson
 */
class FluxReconstructionHyperPoissonModule : public Environment::ModuleRegister<FluxReconstructionHyperPoissonModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionHyperPoisson";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Flux Reconstruction space discretization for Hyperbolized Poisson.";
  }

}; // end FluxReconstructionHyperPoissonModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_FluxReconstructionHyperPoisson_hh

