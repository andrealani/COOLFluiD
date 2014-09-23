#ifndef COOLFluiD_SpectralFDNavierStokes_hh
#define COOLFluiD_SpectralFDNavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Spectral Finite Difference space discretization.
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SpectralFDNavierStokes
 */
class SpectralFDNavierStokesModule : public Environment::ModuleRegister<SpectralFDNavierStokesModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SpectralFDNavierStokes";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the classes specific for the Spectral Finite Difference space discretization for Navier-Stokes.";
  }

}; // end SpectralFDNavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_SpectralFDNavierStokes_hh

