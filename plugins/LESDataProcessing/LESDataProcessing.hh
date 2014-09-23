#ifndef COOLFluiD_Numerics_LESDataProcessing_hh
#define COOLFluiD_Numerics_LESDataProcessing_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    /** 
     * The classes that implement the LESDataProcessing method
     */
    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module LESDataProcessing
 * 
 * @author Willem Deconinck
 */
class LESDataProcessingModule : public Environment::ModuleRegister<LESDataProcessingModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "LESDataProcessing";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the LESDataProcessing method.";
  }

}; // end LESDataProcessingModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESDataProcessing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_LESDataProcessing_hh

