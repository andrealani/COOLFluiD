#ifndef COOLFluiD_Physics_LESvki_hh
#define COOLFluiD_Physics_LESvki_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace LESvki {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module LESvki
 */
class LESvkiModule : public Environment::ModuleRegister<LESvkiModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "LESvki";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a LES flow physical model.";
  }

}; // end LESvkiModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESvki

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_NavierStokes_hh
