#ifndef COOLFluiD_Physics_EntropyMHD_hh
#define COOLFluiD_Physics_EntropyMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement the EntropyMHD MHD model with GLM projection.
  namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module EntropyMHD
 */
class EntropyMHDModule : public Environment::ModuleRegister<EntropyMHDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "EntropyMHD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the EntropyMHD MHD model with GLM divergence cleaning.";
  }

}; // end EntropyMHDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_EntropyMHD_hh
