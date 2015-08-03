#ifndef COOLFluiD_Numerics_FiniteVolumeLTE_hh
#define COOLFluiD_Numerics_FiniteVolumeLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

   namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeLTE
 */
class FiniteVolumeLTEModule : public Environment::ModuleRegister<FiniteVolumeLTEModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeLTE";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the LTE bindings for the FiniteVolume space discretization method.";
  }

}; // end FiniteVolumeLTEModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeLTE_hh

