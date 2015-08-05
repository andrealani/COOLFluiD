#ifndef COOLFluiD_Numerics_FiniteVolumeRadiation_hh
#define COOLFluiD_Numerics_FiniteVolumeRadiation_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeRadiation {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeRadiation
 */
class FiniteVolumeRadiationModule : public Environment::ModuleRegister<FiniteVolumeRadiationModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeRadiation";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module concerns Inductively Coupled Plasma : FiniteVolumeRadiation.";
  }

}; // end FiniteVolumeRadiationModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeRadiation

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeRadiation_hh
