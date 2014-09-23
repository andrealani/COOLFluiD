#ifndef COOLFluiD_Numerics_FiniteVolumeCombustion_hh
#define COOLFluiD_Numerics_FiniteVolumeCombustion_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

   namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeCombustion
 */
class FiniteVolumeCombustionModule : public Environment::ModuleRegister<FiniteVolumeCombustionModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeCombustion";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the bindings for the Combustion model in the FiniteVolume space discretization method.";
  }

}; // end FiniteVolumeCombustionModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeCombustion_hh

