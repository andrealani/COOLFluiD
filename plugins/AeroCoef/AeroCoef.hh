#ifndef COOLFluiD_Numerics_AeroCoef_hh
#define COOLFluiD_Numerics_AeroCoef_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement computation of aerodynamic coefficients
  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module AeroCoef
 */
class AeroCoefModule : public Environment::ModuleRegister<AeroCoefModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "AeroCoef";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements  computation of aerodynamic coefficients.";
  }

}; // end AeroCoefModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_AeroCoef_hh

