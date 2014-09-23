#ifndef COOLFluiD_Numerics_SubSystemCoupler_hh
#define COOLFluiD_Numerics_SubSystemCoupler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement a SubSystem coupler
  namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SubSystemCoupler
 */
class SubSystemCouplerModule : public Environment::ModuleRegister<SubSystemCouplerModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SubSystemCoupler";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a SubSystem coupler.";
  }

}; // end SubSystemCouplerModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_SubSystemCoupler_hh

