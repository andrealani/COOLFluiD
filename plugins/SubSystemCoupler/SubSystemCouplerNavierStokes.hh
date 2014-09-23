#ifndef COOLFluiD_Numerics_SubSystemCouplerNavierStokes_hh
#define COOLFluiD_Numerics_SubSystemCouplerNavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement a SubSystem coupler
  namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SubSystemCouplerNavierStokes
 */
class SubSystemCouplerNavierStokesModule : public Environment::ModuleRegister<SubSystemCouplerNavierStokesModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SubSystemCouplerNavierStokes";
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

}; // end SubSystemCouplerNavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_SubSystemCouplerNavierStokes_hh

