#ifndef COOLFluiD_Physics_NavierStokes_hh
#define COOLFluiD_Physics_NavierStokes_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module NavierStokes
 */
class NavierStokesModule : public Environment::ModuleRegister<NavierStokesModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "NavierStokes";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a NavierStokes flow physical model.";
  }

}; // end NavierStokesModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_NavierStokes_hh
