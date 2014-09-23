#ifndef COOLFluiD_Physics_NitrogenNASA_hh
#define COOLFluiD_Physics_NitrogenNASA_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    /// The classes that implement an interface to the NitrogenNASA library.
    namespace NitrogenNASA {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module NitrogenNASA
 */
class NitrogenNASAModule : public Environment::ModuleRegister<NitrogenNASAModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "NitrogenNASA";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the NitrogenNASA library.";
  }

}; // end NitrogenNASAModule

//////////////////////////////////////////////////////////////////////////////

    }  // namespace NitrogenNASA

  }  // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NitrogenNASA_hh
