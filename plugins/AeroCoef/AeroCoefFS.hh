#ifndef COOLFluiD_Numerics_AeroCoefFS_hh
#define COOLFluiD_Numerics_AeroCoefFS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement computation of aerodynamic coefficients
    /// with bindings to FVM methods
  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module AeroCoefFS
 */
class AeroCoefFSModule : public Environment::ModuleRegister<AeroCoefFSModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "AeroCoefFS";
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

}; // end AeroCoefFSModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_AeroCoefFS_hh

