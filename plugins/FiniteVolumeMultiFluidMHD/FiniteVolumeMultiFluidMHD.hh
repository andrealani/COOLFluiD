#ifndef COOLFluiD_Numerics_FiniteVolumeMultiFluidMHD_hh
#define COOLFluiD_Numerics_FiniteVolumeMultiFluidMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

   namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeMultiFluidMHD
 */
class FiniteVolumeMultiFluidMHDModule : public Environment::ModuleRegister<FiniteVolumeMultiFluidMHDModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeMultiFluidMHD";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements the bindings for the MultiFluidMHD model in the FiniteVolume space discretization method.";
  }

}; // end FiniteVolumeMultiFluidMHDModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeMultiFluidMHD_hh

