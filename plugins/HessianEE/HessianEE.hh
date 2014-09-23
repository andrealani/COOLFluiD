#ifndef COOLFluiD_Numerics_HessianEE_hh
#define COOLFluiD_Numerics_HessianEE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement an Hessian based Error estimator
  namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module HessianEE
 * @author Tiago Quintino
 */
class HessianEEModule : public Environment::ModuleRegister<HessianEEModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "HessianEE";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an Hessian based Error estimator.";
  }

}; // end HessianEEModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_HessianEE_hh

