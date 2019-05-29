#ifndef COOLFluiD_Numerics_FiniteVolumeAdvectionDiffusion_hh
#define COOLFluiD_Numerics_FiniteVolumeAdvectionDiffusion_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    /// Classes concerning AdvectionDiffusion (FiniteVolumeAdvectionDiffusion)
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeAdvectionDiffusion
 */
class FiniteVolumeAdvectionDiffusionModule : public Environment::ModuleRegister<FiniteVolumeAdvectionDiffusionModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeAdvectionDiffusion";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module concerns AdvectionDiffusion Equation solver : FiniteVolumeAdvectionDiffusion.";
  }

}; // end FiniteVolumeAdvectionDiffusionModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeAdvectionDiffusion_hh
