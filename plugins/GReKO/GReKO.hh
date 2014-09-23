#ifndef COOLFluiD_Physics_GReKO_hh
#define COOLFluiD_Physics_GReKO_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a GReKO flow physical model.
  namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module GReKO
 *
 * @author Khalil Bensassi
  */


class GReKOModule : public Environment::ModuleRegister<GReKOModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "GReKO";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a GReKO flow physical model.";
  }

}; // end GReKOModule

//////////////////////////////////////////////////////////////////////////////

    } // end namespace GReKO

  } // end namespace Physics

} // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_GReKO_hh
