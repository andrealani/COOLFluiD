#ifndef COOLFLUID_Numerics_ConcurrentCoupler_ConcurrentCoupler_hh
#define COOLFLUID_Numerics_ConcurrentCoupler_ConcurrentCoupler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
    
    /// The classes that implement a SubSystem coupler
    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module ConcurrentCoupler
 */
class ConcurrentCouplerModule : public Environment::ModuleRegister<ConcurrentCouplerModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ConcurrentCoupler";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a coupler for subsystems running concurrently.";
  }
  
}; // end ConcurrentCouplerMethod

//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_ConcurrentCoupler_ConcurrentCoupler_hh

