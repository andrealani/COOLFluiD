#ifndef COOLFluiD_Physics_ICP_hh
#define COOLFluiD_Physics_ICP_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a ICP flow physical model.
  namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module ICP
 */
class ICPModule : public Environment::ModuleRegister<ICPModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "ICP";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a ICP flow physical model.";
  }

}; // end ICPModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_ICP_hh
