#ifndef COOLFluiD_Post_Process_Module_hh
#define COOLFluiD_Post_Process_Module_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement Grey Media calculation
  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module PostProcess
 */
class PostProcessModule : public Environment::ModuleRegister<PostProcessModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "PostProcess";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface for a post-process routine";
  }

}; // end PostProcessModule

//////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Post_Process_Module_hh
