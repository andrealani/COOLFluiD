#ifndef COOLFluiD_Numerics_SubSystemCouplerHeat_hh
#define COOLFluiD_Numerics_SubSystemCouplerHeat_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement a SubSystem coupler
  namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SubSystemCouplerHeat
 */
class SubSystemCouplerHeatModule : public Environment::ModuleRegister<SubSystemCouplerHeatModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SubSystemCouplerHeat";
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

}; // end SubSystemCouplerHeatModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_SubSystemCouplerHeat_hh

