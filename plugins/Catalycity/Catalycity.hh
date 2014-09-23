#ifndef COOLFluiD_Numerics_Catalycity_hh
#define COOLFluiD_Numerics_Catalycity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement computation of aerodynamic coefficients
  namespace Catalycity {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Catalycity
 */
class CatalycityModule : public Environment::ModuleRegister<CatalycityModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Catalycity";
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

}; // end CatalycityModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Catalycity

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_Catalycity_hh

