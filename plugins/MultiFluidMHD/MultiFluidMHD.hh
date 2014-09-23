#ifndef COOLFluiD_Physics_MultiFluidMHD_hh
#define COOLFluiD_Physics_MultiFluidMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement Ideal MultiFluidMHD flow physical model.
  namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MultiFluidMHD
 */
class MultiFluidMHDModule : public Environment::ModuleRegister<MultiFluidMHDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MultiFluidMHD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a MultiFluidMHD physical model.";
  }

}; // end MultiFluidMHDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_MultiFluidMHD_hh
