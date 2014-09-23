#ifndef COOLFluiD_Physics_MHD_hh
#define COOLFluiD_Physics_MHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement Ideal MHD flow physical model.
  namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MHD
 */
class MHDModule : public Environment::ModuleRegister<MHDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MHD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a Ideal MHD flow physical model.";
  }

}; // end MHDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_MHD_hh
