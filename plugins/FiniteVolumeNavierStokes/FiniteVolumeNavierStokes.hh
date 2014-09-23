#ifndef COOLFluiD_Numerics_FiniteVolumeNavierStokes_hh
#define COOLFluiD_Numerics_FiniteVolumeNavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

   namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeNavierStokes
 */
class FiniteVolumeNavierStokesModule : public Environment::ModuleRegister<FiniteVolumeNavierStokesModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeNavierStokes";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the NavierStokes bindings for the FiniteVolume space discretization method.";
  }

}; // end FiniteVolumeNavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeNavierStokes_hh

