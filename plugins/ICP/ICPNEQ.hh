#ifndef COOLFluiD_Physics_ICPNEQ_hh
#define COOLFluiD_Physics_ICPNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

  namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module ICPNEQ
 */
class ICPNEQModule : public Environment::ModuleRegister<ICPNEQModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ICPNEQ";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a ICPNEQ flow physical model.";
  }

}; // end ICPNEQModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_ICPNEQ_hh
