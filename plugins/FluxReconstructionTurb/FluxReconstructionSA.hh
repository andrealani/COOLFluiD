#ifndef COOLFluiD_FluxReconstructionSA_hh
#define COOLFluiD_FluxReconstructionSA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

   namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionSA
 */
class FluxReconstructionSAModule : public Environment::ModuleRegister<FluxReconstructionSAModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionSA";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the bindings for the Sparlat-Almaras turbulence model in the FluxReconstruction space discretization method.";
  }

}; // end FluxReconstructionSAModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstruction

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FluxReconstructionSA_hh

