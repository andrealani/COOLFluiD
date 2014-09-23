#ifndef COOLFluiD_SpectralFVNavierStokes_hh
#define COOLFluiD_SpectralFVNavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Spectral Finite Volume space discretization.
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SpectralFVNavierStokes
 */
class SpectralFVNavierStokesModule : public Environment::ModuleRegister<SpectralFVNavierStokesModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SpectralFVNavierStokes";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Spectral Finite Volume space discretization for Navier-Stokes.";
  }

}; // end SpectralFVNavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_SpectralFVNavierStokes_hh

