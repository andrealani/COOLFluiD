#ifndef COOLFluiD_FluxReconstructionGammaAlpha_hh
#define COOLFluiD_FluxReconstructionGammaAlpha_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FluxReconstructionGammaAlpha
 *
 *@author Ray Vandenhoeck
 */
class FluxReconstructionGammaAlphaModule :
       public Environment::ModuleRegister<FluxReconstructionGammaAlphaModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FluxReconstructionGammaAlpha";
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

}; // end FluxReconstructionGammaAlphaModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstruction

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_FluxReconstructionGammaAlpha_hh
