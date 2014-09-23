#ifndef COOLFluiD_Numerics_GradComputerFS_hh
#define COOLFluiD_Numerics_GradComputerFS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement computation of aerodynamic coefficients
    /// with bindings to FVM methods
  namespace GradComputer {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module GradComputerFS
 */
class GradComputerFSModule : public Environment::ModuleRegister<GradComputerFSModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "GradComputerFS";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements computation of aerodynamic coefficients with bindings to FS methods.";
  }

}; // end GradComputerFSModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace GradComputer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_GradComputerFS_hh

