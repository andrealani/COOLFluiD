#ifndef COOLFluiD_Numerics_AeroCoefDG_hh
#define COOLFluiD_Numerics_AeroCoefDG_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement computation of aerodynamic coefficients
    /// with bindings to DG methods
  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module AeroCoefDG
 */
class AeroCoefDGModule : public Environment::ModuleRegister<AeroCoefDGModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName() {  return "AeroCoefDG"; }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements computation of aerodynamic coefficients with bindings to DG methods.";
  }

}; // end AeroCoefDGModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_AeroCoefDG_hh

