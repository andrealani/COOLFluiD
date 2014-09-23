#ifndef COOLFluiD_Numerics_LUSGSMethod_hh
#define COOLFluiD_Numerics_LUSGSMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement the (nonlinear) LU-SGS method time stepper.
  namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module LUSGSMethod
 */
class LUSGSMethodModule : public Environment::ModuleRegister<LUSGSMethodModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "LUSGSMethod";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a (nonlinear) LU-SGS method time stepper.";
  }

}; // end LUSGSMethodModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_LUSGSMethod_hh
