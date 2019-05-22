#ifndef COOLFluiD_FluxReconstructionKOmega_hh
#define COOLFluiD_FluxReconstructionKOmega_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionKOmega
 */
class FluxReconstructionKOmegaModule :
       public Environment::ModuleRegister<FluxReconstructionKOmegaModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionKOmega";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the bindings for the K-Omega turbulence model in the FluxReconstruction space discretization method.";
  }

}; // end FluxReconstructionKOmegaModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstruction

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FluxReconstructionKOmega_hh
