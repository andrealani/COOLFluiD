#ifndef COOLFluiD_AutoTemplateLoader_AutoTemplateLoader_hh
#define COOLFluiD_AutoTemplateLoader_AutoTemplateLoader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// This module is still under developement
  namespace AutoTemplateLoader {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module AutoTemplateLoader
 *
 * @author Tiago Quintino
 */
class AutoTemplateLoaderModule : public Environment::ModuleRegister<AutoTemplateLoaderModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "AutoTemplateLoader";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module is still under developement";
  }

}; // end AutoTemplateLoaderModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace AutoTemplateLoader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_AutoTemplateLoader_AutoTemplateLoader_hh

