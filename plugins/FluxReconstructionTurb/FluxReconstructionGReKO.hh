#ifndef COOLFluiD_FluxReconstructionGReKO_hh
#define COOLFluiD_FluxReconstructionGReKO_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionGReKO
 *
 *@author Ray Vandenhoeck
 */
class FluxReconstructionGReKOModule :
       public Environment::ModuleRegister<FluxReconstructionGReKOModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionGReKO";
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

}; // end FluxReconstructionGReKOModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstruction

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_FluxReconstructionGReKO_hh
