#ifndef COOLFluiD_Numerics_FiniteVolumeICP_hh
#define COOLFluiD_Numerics_FiniteVolumeICP_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    /// Classes concerning Inductively Coupled Plasma (FiniteVolumeICP)
    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FiniteVolumeICP
 */
class FiniteVolumeICPModule : public Environment::ModuleRegister<FiniteVolumeICPModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FiniteVolumeICP";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module concerns Inductively Coupled Plasma : FiniteVolumeICP.";
  }

}; // end FiniteVolumeICPModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_FiniteVolumeICP_hh
