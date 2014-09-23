#ifndef COOLFluiD_Numerics_FiniteVolumeGReKO_hh
#define COOLFluiD_Numerics_FiniteVolumeGReKO_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeGReKO
 *
 *@author Khalil Bensassi 
 */
class FiniteVolumeGReKOModule :
       public Environment::ModuleRegister<FiniteVolumeGReKOModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeGReKO";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the bindings for the K-Omega turbulence model in the FiniteVolume space discretization method.";
  }

}; // end FiniteVolumeGReKOModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeGReKO_hh
